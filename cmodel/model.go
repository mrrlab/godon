// Package cmodel provides codon evolution models.
package cmodel

import (
	"math"
	"math/rand"
	"runtime"
	"strconv"
	"sync"

	"github.com/gonum/blas"
	"github.com/gonum/blas/native"
	"github.com/gonum/matrix/mat64"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/tree"
)

// impl provides a type of blas implementation.
var impl native.Implementation

const (
	// If the proportion of site class is less than this number no
	// need to compute probability.
	smallProp = 1e-20
	// minBrLen is the minimal branch length for optimization.
	minBrLen = 1e-9
	// Default value for the maximum branch length.
	defaultMaxBrLen = 100
)

// TreeOptimizable is an extension of optimize.Optimizable which
// includes tree-related methods.
type TreeOptimizable interface {
	optimize.Optimizable
	// SetOptimizeBranchLengths enables branch-length optimization.
	SetOptimizeBranchLengths()
	// GetOptimizeBranchLengths returns true if branch-length
	// optimization is enabled.
	GetOptimizeBranchLengths() bool
	// SetAdaptive enables adaptive MCMC for the TreeOptimizable.
	SetAdaptive(*optimize.AdaptiveSettings)
	// SetMaxBranchLength changes the maximum branch length for
	// the optimization.
	SetMaxBranchLength(float64)
	// SetAggregationMode changes the aggregation mode.
	SetAggregationMode(AggMode)
	// GetTreeString returns tree in a newick format.
	GetTreeString() string
	// Final performs analysis after optimization is complete.
	Final()
	// Summary returns summary of the object for JSON export.
	Summary() interface{}
}

// TreeOptimizableSiteClass is a special case of TreeOptimizable which
// has site-classes.
type TreeOptimizableSiteClass interface {
	TreeOptimizable
	// GetNClass returns number of site classes.
	GetNClass() int
}

// Model is an interface for the model. It provides information about
// number of classes and allows adding parametes.
type Model interface {
	// GetNClass returns number of site classes.
	GetNClass() int
	// addParameters all the parameters of the Model.
	addParameters(optimize.FloatParameterGenerator)
}

// BaseModel stores tree and alignment. Matrices and site classes are
// stored and cached as well.
type BaseModel struct {
	// Model is the model implementation.
	Model

	// data stores codon model input data (i.e. tree, alignment & gcode)
	data *Data

	optBranch bool
	lettersF  [][]int
	lettersA  [][]int
	rshuffle  []int //random shuffle of positions
	//precomtputed aggregation schemas
	schemas []*aggSchema
	qs      [][]*codon.EMatrix
	scale   []float64
	//prop is proportions, in theory it can be different for every site
	//by default it's the same, i.e. prop[0] == prop[1] == ... = prop[ncodons]
	//if it is different special care should be taken in Copy method of the model
	prop       [][]float64
	nclass     int
	parameters optimize.FloatParameters
	as         *optimize.AdaptiveSettings
	maxBrLen   float64

	// aggregation mode
	aggMode AggMode

	// remember computations wee need to perform
	expAllBr bool
	expBr    []bool

	// this is a list of exponentiated matrices
	eQts [][][]float64

	// likelihoods per position
	prunAllPos bool
	prunPos    []bool
	l          []float64
}

// NewBaseModel creates a new base Model.
func NewBaseModel(data *Data, model Model) (bm *BaseModel) {
	f, a := data.cSeqs.Letters()
	nclass := model.GetNClass()
	bm = &BaseModel{
		Model:    model,
		data:     data,
		lettersF: f,
		lettersA: a,
		rshuffle: rand.Perm(data.cSeqs.Length()),
		schemas:  make([]*aggSchema, data.cSeqs.Length()),
		qs:       make([][]*codon.EMatrix, nclass),
		scale:    make([]float64, data.Tree.MaxNodeID()+1),
		expBr:    make([]bool, data.Tree.MaxNodeID()+1),
		prop:     make([][]float64, data.cSeqs.Length()),
		nclass:   nclass,
		l:        make([]float64, data.cSeqs.Length()),
		prunPos:  make([]bool, data.cSeqs.Length()),
	}
	p := make([]float64, nclass)
	for i := range bm.prop {
		bm.prop[i] = p
	}
	for i := 0; i < nclass; i++ {
		bm.qs[i] = make([]*codon.EMatrix, data.Tree.MaxNodeID()+1)
	}
	data.Tree.NodeOrder()
	bm.ReorderAlignment()
	return
}

// Copy creates a copy of BaseModel.
func (m *BaseModel) Copy() (newM *BaseModel) {
	newM = NewBaseModel(m.data.Copy(), m.Model)
	copy(newM.prop[0], m.prop[0])
	newM.as = m.as
	newM.optBranch = m.optBranch
	newM.rshuffle = m.rshuffle
	return
}

// SetAdaptive enables adaptive mode (for adaptive MCMC).
func (m *BaseModel) SetAdaptive(as *optimize.AdaptiveSettings) {
	m.as = as
	m.setupParameters()
}

// GetOptimizeBranchLengths returns true if branch-length
// optimization is enabled.
func (m *BaseModel) GetOptimizeBranchLengths() bool {
	return m.optBranch
}

// SetOptimizeBranchLengths enables branch-length optimization.
func (m *BaseModel) SetOptimizeBranchLengths() {
	m.optBranch = true
	m.setupParameters()
}

// SetMaxBranchLength changes the maximum branch length for
// the optimization.
func (m *BaseModel) SetMaxBranchLength(maxBrLen float64) {
	m.maxBrLen = maxBrLen
	m.setupParameters()
}

// GetFloatParameters returns all the optimization parameters.
func (m *BaseModel) GetFloatParameters() optimize.FloatParameters {
	return m.parameters
}

// setupParameters first delete all the parameters and then adds
// them. This is useful after setting adaptive MCMC mode or other
// changes in the parameters.
func (m *BaseModel) setupParameters() {
	m.parameters = nil
	var fpg optimize.FloatParameterGenerator
	if m.as != nil {
		fpg = m.as.ParameterGenerator
	} else {
		fpg = optimize.BasicFloatParameterGenerator
	}
	m.addBranchParameters(fpg)
	m.Model.addParameters(fpg)
}

// Make branch length parameters adaptive.
func (m *BaseModel) addBranchParameters(fpg optimize.FloatParameterGenerator) {
	if m.maxBrLen == 0 {
		m.maxBrLen = defaultMaxBrLen
	}
	if m.optBranch {

		for _, node := range m.data.Tree.NodeIDArray() {
			if node == nil {
				continue
			}
			nodeID := node.ID
			// Root branch is not optimized
			if node.IsRoot() {
				continue
			}
			par := fpg(&node.BranchLength, "br"+strconv.Itoa(node.ID))
			par.SetOnChange(func() {
				m.expBr[nodeID] = false
			})
			par.SetPriorFunc(optimize.GammaPrior(1, 2, false))
			par.SetMin(minBrLen)
			par.SetMax(m.maxBrLen)
			par.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(par)

		}
	}
}

// SetAggregationMode changes the aggregation mode.
func (m *BaseModel) SetAggregationMode(mode AggMode) {
	m.aggMode = mode
}

// ReorderAlignment reorders codon alignment so order of nodes and
// sequences are the same.  This allows faster access to sequences by
// their index in the array.
func (m *BaseModel) ReorderAlignment() {
	nm2id := make(map[string]int)
	for i, s := range m.data.cSeqs {
		nm2id[s.Name] = i
	}

	if m.data.Tree.NLeaves() != len(m.data.cSeqs) {
		log.Fatal("Tree doesn't match the alignment.")
	}
	newCali := make(codon.Sequences, m.data.Tree.NLeaves())
	for node := range m.data.Tree.Terminals() {
		nodeID, ok := nm2id[node.Name]
		if !ok {
			log.Fatalf("No sequence found for the leaf <%s>.", node.Name)
		}
		newCali[node.LeafID] = m.data.cSeqs[nodeID]
	}

	m.data.cSeqs = newCali
}

// ExpBranch exponentiates a signle branch. This uses eigen decomposed matrices.
func (m *BaseModel) ExpBranch(br int) {
	node := m.data.Tree.NodeIDArray()[br]
	cD := mat64.NewDense(m.data.cFreq.GCode.NCodon, m.data.cFreq.GCode.NCodon, nil)
	for class := range m.qs {
		var oclass int
		for oclass = class - 1; oclass >= 0; oclass-- {
			if m.qs[class][node.ID] == m.qs[oclass][node.ID] {
				m.eQts[class][node.ID] = m.eQts[oclass][node.ID]
				break
			}
		}
		if oclass < 0 {
			Q, err := m.qs[class][node.ID].Exp(cD, node.BranchLength/m.scale[node.ID])
			if err != nil {
				panic("Error exponentiating")
			}
			m.eQts[class][node.ID] = Q.RawMatrix().Data
		}
	}
	m.expBr[br] = true
	m.prunAllPos = false
}

// ExpBranches sxponentiates all branches in the tree.
func (m *BaseModel) ExpBranches() {
	if m.eQts == nil {
		m.eQts = make([][][]float64, len(m.qs))
		for class := range m.qs {
			m.eQts[class] = make([][]float64, m.data.Tree.MaxNodeID()+1)
		}
	} else {
		for class := range m.eQts {
			for nd := range m.eQts[class] {
				m.eQts[class][nd] = nil
			}
		}
	}

	nTasks := len(m.qs) * m.data.Tree.NNodes()

	// expTask is a type storing a task of exponentiating matrices for a
	// class & a node.
	type expTask struct {
		class int
		node  *tree.Node
	}

	tasks := make(chan expTask, nTasks)
	var wg sync.WaitGroup

	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		wg.Add(1)
		go func() {
			cD := mat64.NewDense(m.data.cFreq.GCode.NCodon, m.data.cFreq.GCode.NCodon, nil)
			for s := range tasks {
				Q, err := m.qs[s.class][s.node.ID].Exp(cD, s.node.BranchLength/m.scale[s.node.ID])
				if err != nil {
					panic("error exponentiating matrix")
				}
				m.eQts[s.class][s.node.ID] = Q.RawMatrix().Data
			}
			wg.Done()
		}()
	}

	for class := range m.qs {
		for _, node := range m.data.Tree.NodeIDArray() {
			if node == nil {
				continue
			}
			var oclass int
			for oclass = class - 1; oclass >= 0; oclass-- {
				if m.qs[class][node.ID] == m.qs[oclass][node.ID] {
					defer func(class, oclass, nid int) {
						m.eQts[class][nid] = m.eQts[oclass][nid]
					}(class, oclass, node.ID)
					break
				}
			}
			if oclass < 0 {
				tasks <- expTask{class, node}
			}
			m.expBr[node.ID] = true
		}
	}
	close(tasks)
	wg.Wait()
	m.expAllBr = true
}

// expBranchesIfNeeded performes matrix exponentiation only if it is
// needed. It should be called before the likelihood computations.
func (m *BaseModel) expBranchesIfNeeded() {
	if !m.expAllBr {
		m.ExpBranches()
		m.prunAllPos = false
	} else {
		for _, node := range m.data.Tree.NodeIDArray() {
			if node == nil {
				continue
			}
			if !m.expBr[node.ID] && node != nil {
				m.ExpBranch(node.ID)
				m.prunAllPos = false
			}
		}
	}
}

// GetTreeString returns tree in a newick format.
func (m *BaseModel) GetTreeString() string {
	return m.data.Tree.ClassString()
}

// Likelihood calculates tree likelihood.
func (m *BaseModel) Likelihood() (lnL float64) {
	log.Debugf("x=%v", m.parameters.Values(nil))

	m.expBranchesIfNeeded()

	if len(m.prop) != m.data.cSeqs.Length() {
		panic("incorrect proportion length")
	}

	nPos := m.data.cSeqs.Length()
	nWorkers := runtime.GOMAXPROCS(0)
	done := make(chan struct{}, nWorkers)
	tasks := make(chan int, nPos)

	for i := 0; i < nWorkers; i++ {
		go func() {
			nni := m.data.Tree.MaxNodeID() + 1
			plh := make([][]float64, nni)
			for i := 0; i < nni; i++ {
				plh[i] = make([]float64, m.data.cFreq.GCode.NCodon+1)
			}
			for pos := range tasks {
				if m.prunAllPos && m.prunPos[pos] {
					continue
				}
				res := 0.0
				for class, p := range m.prop[pos] {
					switch {
					case p <= smallProp:
						// if proportion is to small
						continue
					case len(m.lettersF[pos]) == 1:
						// no letters in the current position
						// probability = 1, res += 0
						res += 1 * p
					case m.aggMode == AggFixed && len(m.lettersF[pos]) == 2:
						res += m.fixedSubL(class, pos, plh) * p
					case m.aggMode == AggObserved:
						res += m.observedSubL(class, pos, plh, m.lettersF[pos], m.lettersA[pos]) * p
					case m.aggMode == AggRandom:
						spos := m.rshuffle[pos]
						res += m.observedSubL(class, pos, plh, m.lettersF[spos], m.lettersA[spos]) * p
					case m.aggMode == AggObservedNew:
						schema := m.schemas[pos]
						if schema == nil {
							schema = m.observedStates(m.lettersF[pos], m.lettersA[pos])
							m.schemas[pos] = schema
						}
						res += m.aggSubL(class, pos, plh, schema) * p
					default:
						res += m.fullSubL(class, pos, plh) * p
					}
				}
				m.l[pos] = math.Log(res)
				m.prunPos[pos] = true
			}
			done <- struct{}{}
		}()
	}

	for pos := 0; pos < nPos; pos++ {
		tasks <- pos
	}
	close(tasks)

	// wait for all assignments to finish
	for i := 0; i < nWorkers; i++ {
		<-done
	}

	for i := 0; i < nPos; i++ {
		lnL += m.l[i]
		m.prunPos[i] = true
	}
	m.prunAllPos = true
	if math.IsNaN(lnL) {
		lnL = math.Inf(-1)
	}
	if math.IsInf(lnL, 0) {
		log.Error("Invalid likelihood value, parameters:")
		log.Error(m.parameters.NamesString())
		log.Fatal(m.parameters.ValuesString())
	}
	log.Debugf("L=%v", lnL)
	return
}

// classLikelihood returns an array (slice) of likelihoods for every
// site class, this can be used to perform NEB/BEB (Naive/Bayes
// empirical Bayes) analysis.
func (m *BaseModel) classLikelihoods() (res [][]float64) {
	res = make([][]float64, m.GetNClass())
	nPos := m.data.cSeqs.Length()
	for i := range res {
		res[i] = make([]float64, nPos)
	}

	m.expBranchesIfNeeded()

	nWorkers := runtime.GOMAXPROCS(0)
	done := make(chan struct{}, nWorkers)
	tasks := make(chan int, nPos)

	for i := 0; i < nWorkers; i++ {
		go func() {
			nni := m.data.Tree.MaxNodeID() + 1
			plh := make([][]float64, nni)
			for i := 0; i < nni; i++ {
				plh[i] = make([]float64, m.data.cFreq.GCode.NCodon+1)
			}
			for pos := range tasks {
				for class, p := range m.prop[pos] {
					switch {
					case p <= smallProp:
						// if proportion is too small
						res[class][pos] = 0
						continue
					case len(m.lettersF[pos]) == 1:
						// no letters in the current position
						// probability = 1
						res[class][pos] = 1 * p
					default:
						res[class][pos] = m.fullSubL(class, pos, plh) * p
					}
				}
			}
			done <- struct{}{}
		}()
	}

	for pos := 0; pos < nPos; pos++ {
		tasks <- pos
	}
	close(tasks)

	// wait for all assignments to finish
	for i := 0; i < nWorkers; i++ {
		<-done
	}

	return
}

// NEBPosterior returns array (slice) of posterior probabilities for a
// given classes.
func (m *BaseModel) NEBPosterior(classes map[int]bool) (res []float64) {
	nPos := m.data.cSeqs.Length()
	nClass := m.GetNClass()
	res = make([]float64, nPos)
	// compute class likelihood matrix
	l := m.classLikelihoods()
	for pos := range res {
		numerator := 0.0
		denominator := 0.0
		for cl := 0; cl < nClass; cl++ {
			if classes[cl] {
				numerator += l[cl][pos]
			}
			denominator += l[cl][pos]
		}
		res[pos] = numerator / denominator
	}
	return
}

// PrintPosterior prints results of posterior analysis.
func (m *BaseModel) PrintPosterior(posterior []float64) {
	log.Notice("pos\tcoodn\taa\tp")

	for i, p := range posterior {
		if p > 0.5 {
			bcodon := byte(0)
			for _, seq := range m.data.cSeqs {
				bcodon = seq.Sequence[i]
				if bcodon != codon.NOCODON {
					break
				}
			}

			codon, ok := m.data.cFreq.GCode.NumCodon[bcodon]
			if !ok {
				codon = "NNN"
			}
			aa, ok := m.data.cFreq.GCode.Map[codon]
			if !ok {
				aa = 'X'
			}

			log.Noticef("%v\t%v\t%c\t%0.3f", i+1, codon, aa, p)
		}
	}
}

// fullSubL calculates likelihood for given site class and position.
func (m *BaseModel) fullSubL(class, pos int, plh [][]float64) (res float64) {
	NCodon := m.data.cFreq.GCode.NCodon

	for i := 0; i < m.data.Tree.MaxNodeID()+1; i++ {
		plh[i][0] = math.NaN()
	}

	for node := range m.data.Tree.Terminals() {
		cod := m.data.cSeqs[node.LeafID].Sequence[pos]
		for l := byte(0); l < byte(NCodon); l++ {
			if cod == codon.NOCODON || l == cod {
				plh[node.ID][l] = 1
			} else {
				plh[node.ID][l] = 0
			}
		}
	}

	mul := make([]float64, NCodon)

	for _, node := range m.data.Tree.NodeOrder() {
		for l1 := 0; l1 < NCodon; l1++ {
			plh[node.ID][l1] = 1
		}
		for _, child := range node.ChildNodes() {
			impl.Dgemv(blas.NoTrans, NCodon, NCodon, 1, m.eQts[class][child.ID], NCodon, plh[child.ID], 1, 0, mul, 1)
			for l1 := 0; l1 < NCodon; l1++ {
				plh[node.ID][l1] *= mul[l1]
			}
		}

		if node.IsRoot() {
			for l := 0; l < NCodon; l++ {
				res += m.data.cFreq.Freq[l] * plh[node.ID][l]
			}
			break
		}

	}
	return
}

// Final performs analysis after optimization is complete. This
// implementation does nothing. This method should be implemented.
func (m *BaseModel) Final() {
}

// Summary returns the run summary.
func (m *BaseModel) Summary() interface{} {
	return nil
}

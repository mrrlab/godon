// Package cmodel provides codon evolution models.
package cmodel

import (
	"math"
	"math/rand"
	"runtime"
	"strconv"
	"sync"

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
	// Default value for the maximum branch length.
	defaultMaxBrLen = 100
)

// TreeOptimizable is an extension of optimize.Optimizable which
// includes tree-related methods.
type TreeOptimizable interface {
	optimize.Optimizable
	// SetOptimizeBranchLengths enables branch-length optimization.
	SetOptimizeBranchLengths()
	// SetAdaptive enables adaptive MCMC for the TreeOptimizable.
	SetAdaptive(*optimize.AdaptiveSettings)
	// SetMaxBranchLength changes the maximum branch length for
	// the optimization.
	SetMaxBranchLength(float64)
	// SetAggregationMode changes the aggregation mode.
	SetAggregationMode(AggMode)
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

// Model stores tree and alignment. Matrices and site classes are
// stored and cached as well.
type BaseModel struct {
	// Model is the model implementation.
	Model

	tree      *tree.Tree
	optBranch bool
	cali      codon.CodonSequences
	lettersF  [][]int
	lettersA  [][]int
	rshuffle  []int //random shuffle of positions
	//precomtputed aggregation schemas
	schemas []*aggSchema
	cf      codon.CodonFrequency
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
func NewBaseModel(cali codon.CodonSequences, t *tree.Tree, cf codon.CodonFrequency, model Model) (bm *BaseModel) {
	f, a := cali.Letters()
	nclass := model.GetNClass()
	bm = &BaseModel{
		Model:    model,
		cali:     cali,
		lettersF: f,
		lettersA: a,
		rshuffle: rand.Perm(cali.Length()),
		schemas:  make([]*aggSchema, cali.Length()),
		tree:     t,
		cf:       cf,
		qs:       make([][]*codon.EMatrix, nclass),
		scale:    make([]float64, t.MaxNodeId()+1),
		expBr:    make([]bool, t.MaxNodeId()+1),
		prop:     make([][]float64, cali.Length()),
		nclass:   nclass,
		l:        make([]float64, cali.Length()),
		prunPos:  make([]bool, cali.Length()),
	}
	p := make([]float64, nclass)
	for i := range bm.prop {
		bm.prop[i] = p
	}
	for i := 0; i < nclass; i++ {
		bm.qs[i] = make([]*codon.EMatrix, t.MaxNodeId()+1)
	}
	t.NodeOrder()
	bm.ReorderAlignment()
	return
}

// Copy creates a copy of BaseModel.
func (m *BaseModel) Copy() (newM *BaseModel) {
	newM = NewBaseModel(m.cali, m.tree.Copy(), m.cf, m.Model)
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

		for _, node := range m.tree.NodeIdArray() {
			if node == nil {
				continue
			}
			nodeId := node.Id
			// Root branch is not optimized
			if node.IsRoot() {
				continue
			}
			par := fpg(&node.BranchLength, "br"+strconv.Itoa(node.Id))
			par.SetOnChange(func() {
				m.expBr[nodeId] = false
			})
			par.SetPriorFunc(optimize.GammaPrior(1, 2, false))
			par.SetMin(0)
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
	for i, s := range m.cali {
		nm2id[s.Name] = i
	}

	if m.tree.NLeaves() != len(m.cali) {
		log.Fatal("Tree doesn't match the alignment.")
	}
	newCali := make(codon.CodonSequences, m.tree.NLeaves())
	for node := range m.tree.Terminals() {
		nodeId, ok := nm2id[node.Name]
		if !ok {
			log.Fatalf("No sequence found for the leaf <%s>.", node.Name)
		}
		newCali[node.LeafId] = m.cali[nodeId]
	}

	m.cali = newCali
}

// expTask is a type storing a task of exponentiating matrices for a
// class & a node.
type expTask struct {
	class int
	node  *tree.Node
}

// ExpBranch exponentiates a signle branch. This uses eigen decomposed matrices.
func (m *BaseModel) ExpBranch(br int) {
	node := m.tree.NodeIdArray()[br]
	cD := mat64.NewDense(codon.NCodon, codon.NCodon, nil)
	for class, _ := range m.qs {
		var oclass int
		for oclass = class - 1; oclass >= 0; oclass-- {
			if m.qs[class][node.Id] == m.qs[oclass][node.Id] {
				m.eQts[class][node.Id] = m.eQts[oclass][node.Id]
				break
			}
		}
		if oclass < 0 {
			Q, err := m.qs[class][node.Id].Exp(cD, node.BranchLength/m.scale[node.Id])
			if err != nil {
				panic("Error exponentiating")
			}
			m.eQts[class][node.Id] = Q.RawMatrix().Data
		}
	}
	m.expBr[br] = true
	m.prunAllPos = false
}

// ExpBranches sxponentiates all branches in the tree.
func (m *BaseModel) ExpBranches() {
	if m.eQts == nil {
		m.eQts = make([][][]float64, len(m.qs))
		for class, _ := range m.qs {
			m.eQts[class] = make([][]float64, m.tree.MaxNodeId()+1)
		}
	} else {
		for class, _ := range m.eQts {
			for nd, _ := range m.eQts[class] {
				m.eQts[class][nd] = nil
			}
		}
	}

	nTasks := len(m.qs) * m.tree.NNodes()
	tasks := make(chan expTask, nTasks)
	var wg sync.WaitGroup

	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		wg.Add(1)
		go func() {
			cD := mat64.NewDense(codon.NCodon, codon.NCodon, nil)
			for s := range tasks {
				Q, err := m.qs[s.class][s.node.Id].Exp(cD, s.node.BranchLength/m.scale[s.node.Id])
				if err != nil {
					panic("error exponentiating matrix")
				}
				m.eQts[s.class][s.node.Id] = Q.RawMatrix().Data
			}
			wg.Done()
		}()
	}

	for class, _ := range m.qs {
		for _, node := range m.tree.NodeIdArray() {
			if node == nil {
				continue
			}
			var oclass int
			for oclass = class - 1; oclass >= 0; oclass-- {
				if m.qs[class][node.Id] == m.qs[oclass][node.Id] {
					defer func(class, oclass, nid int) {
						m.eQts[class][nid] = m.eQts[oclass][nid]
					}(class, oclass, node.Id)
					break
				}
			}
			if oclass < 0 {
				tasks <- expTask{class, node}
			}
			m.expBr[node.Id] = true
		}
	}
	close(tasks)
	wg.Wait()
	m.expAllBr = true
}

// Likelihood calculates tree likelihood.
func (m *BaseModel) Likelihood() (lnL float64) {
	log.Debugf("x=%v", m.parameters.Values(nil))
	if !m.expAllBr {
		m.ExpBranches()
		m.prunAllPos = false
	} else {
		for _, node := range m.tree.NodeIdArray() {
			if node == nil {
				continue
			}
			if !m.expBr[node.Id] && node != nil {
				m.ExpBranch(node.Id)
				m.prunAllPos = false
			}
		}
	}

	if len(m.prop) != m.cali.Length() {
		panic("incorrect proportion length")
	}

	nPos := m.cali.Length()
	nWorkers := runtime.GOMAXPROCS(0)
	done := make(chan struct{}, nWorkers)
	tasks := make(chan int, nPos)

	for i := 0; i < nWorkers; i++ {
		go func() {
			nni := m.tree.MaxNodeId() + 1
			plh := make([][]float64, nni)
			for i := 0; i < nni; i++ {
				plh[i] = make([]float64, codon.NCodon+1)
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
					case m.aggMode == AGG_FIXED && len(m.lettersF[pos]) == 2:
						res += m.fixedSubL(class, pos, plh) * p
					case m.aggMode == AGG_OBSERVED:
						res += m.observedSubL(class, pos, plh, m.lettersF[pos], m.lettersA[pos]) * p
					case m.aggMode == AGG_RANDOM:
						spos := m.rshuffle[pos]
						res += m.observedSubL(class, pos, plh, m.lettersF[spos], m.lettersA[spos]) * p
					case m.aggMode == AGG_OBSERVED_NEW:
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
	log.Debugf("L=%v", lnL)
	return
}

// fullSubL calculates likelihood for given site class and position.
func (m *BaseModel) fullSubL(class, pos int, plh [][]float64) (res float64) {
	for i := 0; i < m.tree.MaxNodeId()+1; i++ {
		plh[i][0] = math.NaN()
	}

	for node := range m.tree.Terminals() {
		cod := m.cali[node.LeafId].Sequence[pos]
		for l := byte(0); l < byte(codon.NCodon); l++ {
			if cod == codon.NOCODON || l == cod {
				plh[node.Id][l] = 1
			} else {
				plh[node.Id][l] = 0
			}
		}
	}

	for _, node := range m.tree.NodeOrder() {
		for l1 := 0; l1 < codon.NCodon; l1++ {
			l := 1.0
			for _, child := range node.ChildNodes() {
				// get the row
				q := m.eQts[class][child.Id][l1*codon.NCodon:]
				// get child partial likelhiood
				cplh := plh[child.Id]
				s := impl.Ddot(codon.NCodon, q, 1, cplh, 1)
				l *= s
			}
			plh[node.Id][l1] = l
		}

		if node.IsRoot() {
			for l := 0; l < codon.NCodon; l++ {
				res += m.cf[l] * plh[node.Id][l]
			}
			break
		}

	}
	return
}

// Package cmodel provides codon evolution models.
package cmodel

import (
	"math"
	"runtime"
	"strconv"
	"sync"

	"github.com/gonum/blas/native"
	"github.com/gonum/matrix/mat64"

	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/tree"
)

var impl native.Implementation

type TreeOptimizable interface {
	optimize.Optimizable
	SetOptimizeBranchLengths()
	SetAdaptive(*optimize.AdaptiveSettings)
	// It's a bit confusing name. Here we
	// enable program optimizations.
	SetOptimizations(fixed, all bool)
}

type Model interface {
	GetNClass() int
	addParameters()
	addAdaptiveParameters()
}

// Model stores tree and alignment. Matrices and site classes are stored and cached as well.
type BaseModel struct {
	Model

	tree       *tree.Tree
	optBranch  bool
	cali       CodonSequences
	lettersF   [][]int
	lettersA   [][]int
	cf         CodonFrequency
	qs         [][]*EMatrix
	scale      []float64
	prop       []float64
	nclass     int
	parameters optimize.FloatParameters
	as         *optimize.AdaptiveSettings

	// optimizations
	optFixed bool
	optAll   bool

	// remember computations wee need to perform
	expAllBr bool
	expBr    []bool

	eQts [][][]float64
}

// Creates a new base Model.
func NewBaseModel(cali CodonSequences, t *tree.Tree, cf CodonFrequency, model Model) (bm *BaseModel) {
	f, a := cali.Letters()
	nclass := model.GetNClass()
	bm = &BaseModel{
		Model:    model,
		cali:     cali,
		lettersF: f,
		lettersA: a,
		tree:     t,
		cf:       cf,
		qs:       make([][]*EMatrix, nclass),
		scale:    make([]float64, t.MaxNodeId()+1),
		expBr:    make([]bool, t.MaxNodeId()+1),
		prop:     make([]float64, nclass),
		nclass:   nclass,
	}
	for i := 0; i < nclass; i++ {
		bm.qs[i] = make([]*EMatrix, t.MaxNodeId()+1)
	}
	t.NodeOrder()
	bm.ReorderAlignment()
	return
}

func (m *BaseModel) Copy() (newM *BaseModel) {
	newM = NewBaseModel(m.cali, m.tree.Copy(), m.cf, m.Model)
	copy(newM.prop, m.prop)
	newM.as = m.as
	newM.optBranch = m.optBranch
	return
}

func (m *BaseModel) SetAdaptive(as *optimize.AdaptiveSettings) {
	m.as = as
	m.setupParameters()
}

func (m *BaseModel) SetOptimizeBranchLengths() {
	m.optBranch = true
	m.setupParameters()
}

func (m *BaseModel) GetFloatParameters() optimize.FloatParameters {
	return m.parameters
}

func (m *BaseModel) setupParameters() {
	m.parameters = nil
	m.addBranchParameters()
	if m.as != nil {
		m.Model.addAdaptiveParameters()
	} else {
		m.Model.addParameters()
	}
}

// Make branch length parameters adaptive.
func (m *BaseModel) addBranchParameters() {
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
			if m.as == nil {
				par := optimize.NewBasicFloatParameter(&node.BranchLength, "br"+strconv.Itoa(node.Id))
				par.OnChange = func() {
					m.expBr[nodeId] = false
				}
				par.PriorFunc = optimize.GammaPrior(1, 2, false)
				par.Min = 0
				par.ProposalFunc = optimize.NormalProposal(0.01)
				m.parameters.Append(par)
			} else {
				par := optimize.NewAdaptiveParameter(&node.BranchLength, "br"+strconv.Itoa(node.Id), m.as)
				par.OnChange = func() {
					m.expBr[nodeId] = false
				}
				par.PriorFunc = optimize.GammaPrior(1, 2, false)
				par.Min = 0
				par.ProposalFunc = optimize.NormalProposal(0.01)
				m.parameters.Append(par)
			}

		}
	}
}

func (m *BaseModel) SetOptimizations(fixed, all bool) {
	m.optFixed = fixed
	m.optAll = all
}

// Reorder codon alignment so order of nodes and sequences are the same.
// This allows faster access to sequences by their index in the array.
func (m *BaseModel) ReorderAlignment() {
	nm2id := make(map[string]int)
	for i, s := range m.cali {
		nm2id[s.Name] = i
	}

	newCali := make(CodonSequences, m.tree.NLeaves())
	for node := range m.tree.Terminals() {
		newCali[node.LeafId] = m.cali[nm2id[node.Name]]
	}

	m.cali = newCali
}

type expTask struct {
	class int
	node  *tree.Node
}

// Exponentiate a signle branch. This uses eigen decomposed matrices.
func (m *BaseModel) ExpBranch(br int) {
	node := m.tree.NodeIdArray()[br]
	cD := mat64.NewDense(nCodon, nCodon, nil)
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
}

// Exponentiate all branches in the tree.
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
			cD := mat64.NewDense(nCodon, nCodon, nil)
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

// Calculate tree likelihood.
func (m *BaseModel) Likelihood() (lnL float64) {
	if !m.expAllBr {
		m.ExpBranches()
	} else {
		for _, node := range m.tree.NodeIdArray() {
			if node == nil {
				continue
			}
			if !m.expBr[node.Id] && node != nil {
				m.ExpBranch(node.Id)
			}
		}
	}

	if len(m.prop) != len(m.qs) {
		panic("incorrect proportion length")
	}

	nTasks := len(m.cali[0].Sequence)
	results := make(chan float64, nTasks)
	tasks := make(chan int, nTasks)

	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		go func() {
			nni := m.tree.MaxNodeId() + 1
			plh := make([][]float64, nni)
			for i := 0; i < nni; i++ {
				plh[i] = make([]float64, nCodon+1)
			}
			for pos := range tasks {
				res := 0.0
				for class, p := range m.prop {
					if m.optFixed && len(m.lettersF[pos]) == 2 {
						res += m.fixedSubL(class, pos, plh) * p
					} else {
						if m.optAll {
							res += m.observedSubL(class, pos, plh) * p
						} else {
							res += m.fullSubL(class, pos, plh) * p
						}
					}
				}
				results <- math.Log(res)
			}
		}()
	}

	for pos := 0; pos < nTasks; pos++ {
		tasks <- pos
	}
	close(tasks)

	for i := 0; i < nTasks; i++ {
		dlnL := <-results
		lnL += dlnL
	}
	if math.IsNaN(lnL) {
		lnL = math.Inf(-1)
	}
	return
}

// fullSubL calculates likelihood for given site class and position.
func (m *BaseModel) fullSubL(class, pos int, plh [][]float64) (res float64) {
	for i := 0; i < m.tree.MaxNodeId()+1; i++ {
		plh[i][0] = math.NaN()
	}

	for node := range m.tree.Terminals() {
		codon := m.cali[node.LeafId].Sequence[pos]
		for l := byte(0); l < byte(nCodon); l++ {
			if codon == NOCODON || l == codon {
				plh[node.Id][l] = 1
			} else {
				plh[node.Id][l] = 0
			}
		}
	}

	for _, node := range m.tree.NodeOrder() {
		for l1 := 0; l1 < nCodon; l1++ {
			l := 1.0
			for _, child := range node.ChildNodes() {
				// get the row
				q := m.eQts[class][child.Id][l1*nCodon:]
				// get child partial likelhiood
				cplh := plh[child.Id]
				s := impl.Ddot(nCodon, q, 1, cplh, 1)
				l *= s
			}
			plh[node.Id][l1] = l
		}

		if node.IsRoot() {
			for l := 0; l < nCodon; l++ {
				res += m.cf[l] * plh[node.Id][l]
			}
			break
		}

	}
	return
}

// observedSubL calculates likelihood for given site class and position
// taking into account only visible states.
func (m *BaseModel) observedSubL(class, pos int, plh [][]float64) (res float64) {
	lettersF := m.lettersF[pos]
	lettersA := m.lettersA[pos]
	fabs := 0.0
	for _, l := range lettersA {
		fabs += m.cf[l]
	}

	for node := range m.tree.Terminals() {
		codon := m.cali[node.LeafId].Sequence[pos]
		for _, l := range lettersF {
			if l == int(m.cali[node.LeafId].Sequence[pos]) || codon == NOCODON {
				plh[node.Id][l] = 1
			} else {
				plh[node.Id][l] = 0
			}
		}
		if codon == NOCODON {
			plh[node.Id][nCodon] = 1
		}
	}

	for _, node := range m.tree.NodeOrder() {
		for _, l1 := range lettersF {
			l := 1.0
			for _, child := range node.ChildNodes() {
				// get child partial likelhiood
				cplh := plh[child.Id]
				s := 0.0
				if l1 != nCodon {
					// get the row
					q := m.eQts[class][child.Id][l1*nCodon:]

					for _, l2 := range lettersF {
						//s += q.Get(l1, l2) * plh[child.Id][l2]
						if l2 != nCodon {
							s += q[l2] * cplh[l2]
						} else {
							pia := 0.0
							for _, l2 := range lettersA {
								pia += q[l2]
							}
							s += pia * cplh[l2]
						}
					}
				} else {

					paa := 1.0
					for _, l2 := range lettersF {
						pai := 0.0
						if l2 != nCodon {
							for _, l1 := range lettersA {
								pai += m.cf[l1] * m.eQts[class][child.Id][l1*nCodon+l2]
							}
							pai /= fabs

							paa -= pai
							s += pai * cplh[l2]
						} else {
							s += paa * cplh[l2]
						}
					}
				}
				l *= s
			}
			plh[node.Id][l1] = l
		}

		if node.IsRoot() {
			for _, l := range lettersF {
				if l != nCodon {

					res += m.cf[l] * plh[node.Id][l]
				} else {
					res += fabs * plh[node.Id][l]
				}
			}
			break
		}

	}
	return
}

// fixedSubL calculates likelihood for given site class and position
// if the site is fixed.
func (m *BaseModel) fixedSubL(class, pos int, plh [][]float64) (res float64) {
	for i := 0; i < m.tree.MaxNodeId()+1; i++ {
		plh[i][0] = math.NaN()
	}

	l := int(m.cali[0].Sequence[pos])

	for node := range m.tree.Terminals() {
		plh[node.Id][0] = 1
		plh[node.Id][1] = 0
	}

	for _, node := range m.tree.NodeOrder() {
		plh[node.Id][0] = 1
		plh[node.Id][1] = 1
		for _, child := range node.ChildNodes() {
			p00 := m.eQts[class][child.Id][l*nCodon+l]
			p01 := 1 - p00
			p10 := 0.0
			for l1 := 0; l1 < nCodon; l1++ {
				if l != l1 {
					p10 += m.cf[l1] * m.eQts[class][child.Id][l1*nCodon+l]
				}
			}
			p10 /= (1 - m.cf[l])
			p11 := 1 - p10

			cplh := plh[child.Id]
			plh[node.Id][0] *= p00*cplh[0] + p01*cplh[1]
			plh[node.Id][1] *= p10*cplh[0] + p11*cplh[1]
		}

		if node.IsRoot() {
			res = m.cf[l]*plh[node.Id][0] + (1-m.cf[l])*plh[node.Id][1]
			break
		}

	}
	return
}

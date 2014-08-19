package main

import (
	"fmt"
	"math"
	"runtime"
	"sync"

	"github.com/skelterjohn/go.matrix"

	"bitbucket.com/Davydov/golh/tree"
)

type Optimizable interface {
	SetDefaults()
	GetNumberOfParameters() int
	GetParameterName(i int) string
	GetParameter(i int) float64
	SetParameter(i int, val float64)
	SetOptBranch(optBranch bool)
	Likelihood() float64
}

type Model struct {
	tree      *tree.Tree
	cali      CodonSequences
	cf        CodonFrequency
	qs        [][]*EMatrix
	scale     []float64
	prop      []float64
	nclass    int
	optBranch bool

	eQts [][][]float64
}

func NewModel(cali CodonSequences, t *tree.Tree, cf CodonFrequency, nclass int) (m *Model) {
	m = &Model{cali: cali,
		tree:   t,
		cf:     cf,
		qs:     make([][]*EMatrix, nclass),
		scale:  make([]float64, t.NNodes()),
		prop:   make([]float64, nclass),
		nclass: nclass,
	}
	for i := 0; i < nclass; i++ {
		m.qs[i] = make([]*EMatrix, t.NNodes())
	}
	t.NodeOrder()
	m.ReorderAlignment()

	return
}

func (m *Model) GetNumberOfParameters() int {
	// root branch is not considered
	if m.optBranch {
		return m.tree.NNodes() - 1
	}
	return 0
}

func (m *Model) GetParameterName(i int) string {
	if !m.optBranch {
		panic("branch length are not optimizable")
	}
	return fmt.Sprintf("br%d", m.tree.Nodes()[i+1].Id)
}

func (m *Model) GetParameter(i int) float64 {
	if !m.optBranch {
		panic("branch length are not optimizable")
	}
	return m.tree.Nodes()[i+1].BranchLength
}

func (m *Model) SetParameter(i int, val float64) {
	if !m.optBranch {
		panic("branch length are not optimizable")
	}
	br := i + 1
	m.tree.Nodes()[br].BranchLength = math.Abs(val)
	m.ExpBranch(br)
}

func (m *Model) SetOptBranch(optBranch bool) {
	m.optBranch = optBranch
}

func (m *Model) ReorderAlignment() {
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

func (m *Model) ExpBranch(br int) {
	node := m.tree.Nodes()[br]
	cD := matrix.Zeros(nCodon, nCodon)
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
			m.eQts[class][node.Id] = Q.Array()
		}
	}
}

func (m *Model) ExpBranches() {
	if m.eQts == nil {
		m.eQts = make([][][]float64, len(m.qs))
		for class, _ := range m.qs {
			m.eQts[class] = make([][]float64, m.tree.NNodes())
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
			cD := matrix.Zeros(nCodon, nCodon)
			for s := range tasks {
				Q, err := m.qs[s.class][s.node.Id].Exp(cD, s.node.BranchLength/m.scale[s.node.Id])
				if err != nil {
					panic("error exponentiating matrix")
				}
				m.eQts[s.class][s.node.Id] = Q.Array()
			}
			wg.Done()
		}()
	}

	for class, _ := range m.qs {
		for _, node := range m.tree.Nodes() {
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
		}
	}
	close(tasks)
	wg.Wait()
}

func (m *Model) Likelihood() (lnL float64) {
	if len(m.prop) != len(m.qs) {
		panic("incorrect proportion length")
	}

	nTasks := len(m.cali[0].Sequence)
	results := make(chan float64, nTasks)
	tasks := make(chan int, nTasks)

	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		go func() {
			plh := make([][]float64, m.tree.NNodes())
			for i := 0; i < m.tree.NNodes(); i++ {
				plh[i] = make([]float64, nCodon)
			}
			for pos := range tasks {
				res := 0.0
				for class, p := range m.prop {
					res += m.subL(class, pos, plh) * p
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
	return
}

func (m *Model) subL(class, pos int, plh [][]float64) (res float64) {
	for i := 0; i < m.tree.NNodes(); i++ {
		plh[i][0] = math.NaN()
	}

	for node := range m.tree.Terminals() {
		for l := byte(0); l < byte(nCodon); l++ {
			if l == m.cali[node.LeafId].Sequence[pos] {
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
				s := 0.0
				for l2 := 0; l2 < nCodon; l2++ {
					//s += q.Get(l1, l2) * plh[child.Id][l2]
					s += q[l2] * cplh[l2]
				}
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

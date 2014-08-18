package main

import (
	"math"
	"runtime"
	"sync"

	"github.com/skelterjohn/go.matrix"

	"bitbucket.com/Davydov/golh/tree"
)

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
			var err error
			m.eQts[class][node.Id], err = m.qs[class][node.Id].Exp(cD, node.BranchLength/m.scale[node.Id])
			if err != nil {
				panic("Error exponentiating")
			}
		}
	}
}

func (m *Model) ExpBranches() {
	if m.eQts == nil {
		m.eQts = make([][]*matrix.DenseMatrix, len(m.qs))
		for class, _ := range m.qs {
			m.eQts[class] = make([]*matrix.DenseMatrix, m.tree.NNodes())
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
			var err error
			cD := matrix.Zeros(nCodon, nCodon)
			for s := range tasks {
				m.eQts[s.class][s.node.Id], err = m.qs[s.class][s.node.Id].Exp(cD, s.node.BranchLength/m.scale[s.node.Id])
				if err != nil {
					panic("error exponentiating matrix")
				}
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
			if l == m.cali[m.nm2id[node.Name]].Sequence[pos] {
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
				s := 0.0
				for l2 := 0; l2 < nCodon; l2++ {
					s += m.eQts[class][child.Id].Get(l1, l2) * plh[child.Id][l2]
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

package main

import (
	"math"
	"sync"

	"github.com/skelterjohn/go.matrix"

	"bitbucket.com/Davydov/golh/tree"
)

type expTask struct {
	class int
	node  *tree.Node
}

func ExpBranch(t *tree.Tree, Qs [][]*EMatrix, scale []float64) (eQts [][]*matrix.DenseMatrix) {
	eQts = make([][]*matrix.DenseMatrix, len(Qs))

	nTasks := len(Qs) * t.NNodes()
	tasks := make(chan expTask, nTasks)
	var wg sync.WaitGroup

	for i := 0; i < mxprc; i++ {
		wg.Add(1)
		go func() {
			var err error
			cD := matrix.Zeros(nCodon, nCodon)
			for s := range tasks {
				eQts[s.class][s.node.Id], err = Qs[s.class][s.node.Id].Exp(cD, s.node.BranchLength/scale[s.node.Id])
				if err != nil {
					panic("error exponentiating matrix")
				}
			}
			wg.Done()
		}()
	}

	for class, _ := range Qs {
		eQts[class] = make([]*matrix.DenseMatrix, t.NNodes())
		for _, node := range t.Nodes() {
			var oclass int
			for oclass = class - 1; oclass >= 0; oclass-- {
				if Qs[class][node.Id] == Qs[oclass][node.Id] {
					defer func(class, oclass, nid int) {
						eQts[class][nid] = eQts[oclass][nid]
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

	return
}

func L(ali CodonSequences, t *tree.Tree, prop []float64, scale []float64, Qs [][]*EMatrix, cf CodonFrequency) (lnL float64) {
	if len(prop) != len(Qs) {
		panic("incorrect proportion length")
	}

	eQts := ExpBranch(t, Qs, scale)

	no := t.NodeOrder()

	nTasks := len(ali[0].Sequence)
	results := make(chan float64, nTasks)
	tasks := make(chan int, nTasks)

	for i := 0; i < mxprc; i++ {
		go func() {
			plh := make([][]float64, t.NNodes())
			for i := 0; i < t.NNodes(); i++ {
				plh[i] = make([]float64, nCodon)
			}
			for pos := range tasks {
				res := float64(0)
				for i, p := range prop {
					res += subL(ali, t, no, eQts[i], cf, pos, plh) * p
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

func subL(ali CodonSequences, t *tree.Tree, no []*tree.Node, eQts []*matrix.DenseMatrix, cf CodonFrequency, pos int, plh [][]float64) (res float64) {

	for i := 0; i < t.NNodes(); i++ {
		plh[i][0] = math.NaN()
	}

	for node := range t.Terminals() {
		for l := byte(0); l < byte(nCodon); l++ {
			if l == ali[nm2id[node.Name]].Sequence[pos] {
				plh[node.Id][l] = 1
			} else {
				plh[node.Id][l] = 0
			}
		}
	}

	for _, node := range no {
		for l1 := 0; l1 < nCodon; l1++ {
			l := 1.0
			for _, child := range node.ChildNodes() {
				s := 0.0
				for l2 := 0; l2 < nCodon; l2++ {
					s += eQts[child.Id].Get(l1, l2) * plh[child.Id][l2]
				}
				l *= s
			}
			plh[node.Id][l1] = l
		}

		if node.IsRoot() {
			for l := 0; l < nCodon; l++ {
				res += cf[l] * plh[node.Id][l]
			}
			break
		}

	}
	return
}

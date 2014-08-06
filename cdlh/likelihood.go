package main

import (
	"math"
	"runtime"

	"github.com/skelterjohn/go.matrix"

	"bitbucket.com/Davydov/golh/tree"
)

func L(ali CodonSequences, t *tree.Tree, prop []float64, Qs [][]*EMatrix, cf CodonFrequency) (lnL float64) {
	if len(prop) != len(Qs) {
		panic("incorrect proportion length")
	}

	cD := matrix.Zeros(nCodon, nCodon)
	eQts := make([][]*matrix.DenseMatrix, len(prop))

	for i, _ := range prop {
		eQts[i] = make([]*matrix.DenseMatrix, t.NNodes())
		for node := range t.Nodes() {
			var err error
			eQts[i][node.Id], err = Qs[i][node.Id].Exp(cD, node.BranchLength)
			if err != nil {
				panic("Error exponentiating matrix")
			}
		}
	}

	mxprc := runtime.GOMAXPROCS(0)
	plhch := make(chan [][]float64, mxprc)
	for i := 0; i < mxprc; i++ {
		plhch <- nil
	}

	ch := make(chan float64, len(ali[0].Sequence))
	for pos, _ := range ali[0].Sequence {
		go func(pos int) {
			res := float64(0)
			for i, p := range prop {
				res += subL(ali, t, eQts[i], cf, pos, plhch) * p
			}
			ch <- math.Log(res)
		}(pos)
	}

	for _, _ = range ali[0].Sequence {
		dlnL := <-ch
		lnL += dlnL
	}
	close(ch)
	return
}

func subL(ali CodonSequences, t *tree.Tree, eQts []*matrix.DenseMatrix, cf CodonFrequency, pos int, plhch chan [][]float64) float64 {
	res := 0.0
	plh := <-plhch
	nNodes := t.NNodes()
	if plh == nil {
		plh = make([][]float64, nNodes)
		for i := 0; i < nNodes; i++ {
			plh[i] = make([]float64, nCodon)
		}
	}

	for i := 0; i < nNodes; i++ {
		plh[i][0] = math.NaN()
	}

	nodes := make(chan *tree.Tree, len(ali))
	for node := range t.Terminals() {
		for l := byte(0); l < byte(nCodon); l++ {
			if l == ali[nm2id[node.Name]].Sequence[pos] {
				plh[node.Id][l] = 1
			} else {
				plh[node.Id][l] = 0
			}
		}
		nodes <- node.Parent
	}

NodeLoop:
	for node := range nodes {
		for child := range node.ChildNodes() {
			if math.IsNaN(plh[child.Id][0]) {
				nodes <- node
				continue NodeLoop
			}
		}
		if !math.IsNaN(plh[node.Id][0]) {
			continue NodeLoop
		}
		for l1 := 0; l1 < nCodon; l1++ {
			l := 1.0
			for child := range node.ChildNodes() {
				s := 0.0
				for l2 := 0; l2 < nCodon; l2++ {
					s += eQts[child.Id].Get(l1, l2) * plh[child.Id][l2]
				}
				l *= s
			}
			plh[node.Id][l1] = l
		}
		nodes <- node.Parent
		if node.IsRoot() {
			for l := 0; l < nCodon; l++ {
				res += cf[l] * plh[node.Id][l]
			}
			break NodeLoop
		}

	}
	close(nodes)
	plhch <- plh
	return res
}

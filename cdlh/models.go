package main

import (
	"fmt"

	"bitbucket.com/Davydov/golh/tree"
)

func M0(cali CodonSequences, t *tree.Tree, cf CodonFrequency, kappa, omega float64) float64 {
	Q, s := createTransitionMatrix(cf, kappa, omega)
	Qs := make([][]*EMatrix, 1)
	Qs[0] = make([]*EMatrix, t.NSubNodes())
	scale := make([]float64, t.NSubNodes())
	em := NewEMatrix(Q)
	err := em.Eigen()
	if err != nil {
		panic(fmt.Sprintf("error finding eigen: %v", err))
	}
	for i := 0; i < len(Qs[0]); i++ {
		Qs[0][i] = em
		scale[i] = s
	}
	return L(cali, t, []float64{1}, scale, Qs, cf)
}

func H1(cali CodonSequences, t *tree.Tree, cf CodonFrequency, kappa float64, omega0, omega2 float64, p0, p1, p2a, p2b float64) float64 {
	//fmt.Printf("kappa=%f, omega0=%f, omega2=%f, p=[%f, %f, %f, %f]\n", kappa, omega0, omega2, p0, p1, p2a, p2b)
	Q0, s0 := createTransitionMatrix(cf, kappa, omega0)
	Q1, s1 := createTransitionMatrix(cf, kappa, 1)
	Q2, s2 := createTransitionMatrix(cf, kappa, omega2)
	em0 := NewEMatrix(Q0)
	em1 := NewEMatrix(Q1)
	em2 := NewEMatrix(Q2)
	err1 := em0.Eigen()
	err2 := em1.Eigen()
	err3 := em2.Eigen()
	if err1 != nil || err2 != nil || err3 != nil {
		panic(fmt.Sprintf("error finding eigen: %v, %v, %v", err1, err2, err3))
	}
	Qs := make([][]*EMatrix, 4)
	scale := make([]float64, t.NSubNodes())
	for node := range t.Nodes() {
		if node.Class == 0 {
			scale[node.Id] = (p0+p2a)*s0 + (p1+p2b)*s1
		} else {
			scale[node.Id] = p0*s0 + p1*s1 + (p2a+p2b)*s2
		}
	}
	for i := 0; i < len(Qs); i++ {
		Qs[i] = make([]*EMatrix, t.NSubNodes())
		for node := range t.Nodes() {
			switch i {
			case 0:
				Qs[i][node.Id] = em0
			case 1:
				Qs[i][node.Id] = em1
			case 2:
				if node.Class == 0 {
					Qs[i][node.Id] = em0
				} else {
					Qs[i][node.Id] = em2
				}
			case 3:
				if node.Class == 0 {
					Qs[i][node.Id] = em1
				} else {
					Qs[i][node.Id] = em2
				}
			}
		}
	}
	return L(cali, t, []float64{p0, p1, p2a, p2b}, scale, Qs, cf)
}

package main

import (
	"fmt"

	"bitbucket.com/Davydov/golh/tree"
)

type Model interface {
	Likelihood()
}

type M0 struct {
	q *EMatrix
	qs [][]*EMatrix
	cf CodonFrequency
	cali CodonSequences
	tree *tree.Tree
	omega, kappa float64
	scale []float64
	prop []float64
}

func NewM0(cali CodonSequences, t *tree.Tree, cf CodonFrequency) (m0 *M0) {
	m0 = &M0{cali: cali,
		tree:  t,
		cf: cf,
		qs: make([][]*EMatrix, 1),
		scale: make([]float64, t.NNodes()),
		prop: []float64{1},
		q: &EMatrix{},
	}
	m0.qs[0] = make([]*EMatrix, t.NNodes())
	return

}

func (m0 *M0) SetParameters(kappa, omega float64) {
	m0.kappa = kappa
	m0.omega = omega
	m0.UpdateMatrices()
}

func (m0 *M0) UpdateMatrices() {
	Q, s := createTransitionMatrix(m0.cf, m0.kappa, m0.omega, m0.q.Q)
	m0.q.Set(Q)

	err := m0.q.Eigen()
	if err != nil {
		panic(fmt.Sprintf("error finding eigen: %v", err))
	}
	for i := 0; i < len(m0.qs[0]); i++ {
		m0.qs[0][i] = m0.q
		m0.scale[i] = s
	}
}

func (m0 *M0) Likelihood() float64 {
	return L(m0.cali, m0.tree, m0.prop, m0.scale, m0.qs, m0.cf)
}

func H1(cali CodonSequences, t *tree.Tree, cf CodonFrequency, kappa float64, omega0, omega2 float64, p0, p1, p2a, p2b float64) float64 {
	//fmt.Printf("kappa=%f, omega0=%f, omega2=%f, p=[%f, %f, %f, %f]\n", kappa, omega0, omega2, p0, p1, p2a, p2b)
	Q0, s0 := createTransitionMatrix(cf, kappa, omega0, nil)
	Q1, s1 := createTransitionMatrix(cf, kappa, 1, nil)
	Q2, s2 := createTransitionMatrix(cf, kappa, omega2, nil)
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
	scale := make([]float64, t.NNodes())
	for _, node := range t.Nodes() {
		if node.Class == 0 {
			scale[node.Id] = (p0+p2a)*s0 + (p1+p2b)*s1
		} else {
			scale[node.Id] = p0*s0 + p1*s1 + (p2a+p2b)*s2
		}
	}
	for i := 0; i < len(Qs); i++ {
		Qs[i] = make([]*EMatrix, t.NNodes())
		for _, node := range t.Nodes() {
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

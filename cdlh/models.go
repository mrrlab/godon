package main

import (
	"fmt"

	"github.com/skelterjohn/go.matrix"

	"bitbucket.com/Davydov/golh/tree"
)

type Model struct {
	tree  *tree.Tree
	cali  CodonSequences
	cf    CodonFrequency
	qs    [][]*EMatrix
	scale []float64
	prop  []float64

	eQts [][]*matrix.DenseMatrix
}

type BranchData struct {
}

type M0 struct {
	Model
	q            *EMatrix
	omega, kappa float64
}

func NewM0(cali CodonSequences, t *tree.Tree, cf CodonFrequency) (m *M0) {
	m = &M0{
		Model: Model{cali: cali,
			tree:  t,
			cf:    cf,
			qs:    make([][]*EMatrix, 1),
			scale: make([]float64, t.NNodes()),
			prop:  []float64{1},
		},
		q: &EMatrix{},
	}
	m.qs[0] = make([]*EMatrix, t.NNodes())
	t.NodeOrder()
	return
}

func (m *M0) SetParameters(kappa, omega float64) {
	m.kappa = kappa
	m.omega = omega
	m.UpdateMatrices()
}

func (m *M0) UpdateMatrices() {
	Q, s := createTransitionMatrix(m.cf, m.kappa, m.omega, m.q.Q)
	m.q.Set(Q)

	err := m.q.Eigen()
	if err != nil {
		panic(fmt.Sprintf("error finding eigen: %v", err))
	}
	for i := 0; i < len(m.qs[0]); i++ {
		m.qs[0][i] = m.q
		m.scale[i] = s
	}
}

type H1 struct {
	Model
	q0, q1, q2     *EMatrix
	kappa          float64
	omega0, omega2 float64
}

func NewH1(cali CodonSequences, t *tree.Tree, cf CodonFrequency) (m *H1) {
	m = &H1{
		Model: Model{cali: cali,
			tree:  t,
			cf:    cf,
			qs:    make([][]*EMatrix, 4),
			scale: make([]float64, t.NNodes()),
			prop:  make([]float64, 4),
		},
		q0: &EMatrix{},
		q1: &EMatrix{},
		q2: &EMatrix{},
	}
	for i := 0; i < len(m.qs); i++ {
		m.qs[i] = make([]*EMatrix, t.NNodes())
	}
	t.NodeOrder()
	return

}

func (m *H1) SetParameters(kappa float64, omega0, omega2 float64, p0, p1, p2a, p2b float64) {
	m.kappa = kappa
	m.omega0 = omega0
	m.omega2 = omega2
	m.prop[0] = p0
	m.prop[1] = p1
	m.prop[2] = p2a
	m.prop[3] = p2b
	m.UpdateMatrices()
}

func (m *H1) UpdateMatrices() {
	Q0, s0 := createTransitionMatrix(m.cf, m.kappa, m.omega0, m.q0.Q)
	m.q0.Set(Q0)
	Q1, s1 := createTransitionMatrix(m.cf, m.kappa, 1, m.q1.Q)
	m.q1.Set(Q1)
	Q2, s2 := createTransitionMatrix(m.cf, m.kappa, m.omega2, m.q2.Q)
	m.q2.Set(Q2)

	err0 := m.q0.Eigen()
	err1 := m.q1.Eigen()
	err2 := m.q2.Eigen()

	if err0 != nil || err1 != nil || err2 != nil {
		panic(fmt.Sprintf("error finding eigen: %v, %v, %v", err0, err1, err2))
	}

	for _, node := range m.tree.Nodes() {
		if node.Class == 0 {
			m.scale[node.Id] = (m.prop[0]+m.prop[2])*s0 + (m.prop[1]+m.prop[3])*s1
		} else {
			m.scale[node.Id] = m.prop[0]*s0 + m.prop[1]*s1 + (m.prop[2]+m.prop[3])*s2
		}
	}

	for i := 0; i < len(m.qs); i++ {
		for _, node := range m.tree.Nodes() {
			switch i {
			case 0:
				m.qs[i][node.Id] = m.q0
			case 1:
				m.qs[i][node.Id] = m.q1
			case 2:
				if node.Class == 0 {
					m.qs[i][node.Id] = m.q0
				} else {
					m.qs[i][node.Id] = m.q2
				}
			case 3:
				if node.Class == 0 {
					m.qs[i][node.Id] = m.q1
				} else {
					m.qs[i][node.Id] = m.q2
				}
			}
		}
	}
}

package main

import (
	"fmt"

	"github.com/skelterjohn/go.matrix"

	"bitbucket.com/Davydov/golh/tree"
)

type Model struct {
	tree   *tree.Tree
	cali   CodonSequences
	cf     CodonFrequency
	qs     [][]*EMatrix
	scale  []float64
	prop   []float64
	nm2id  map[string]int
	nclass int

	eQts [][]*matrix.DenseMatrix
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

	m.nm2id = make(map[string]int)
	for i, s := range m.cali {
		m.nm2id[s.Name] = i
	}

	return
}

type M0 struct {
	*Model
	q            *EMatrix
	omega, kappa float64
}

func NewM0(cali CodonSequences, t *tree.Tree, cf CodonFrequency) (m *M0) {
	m = &M0{
		Model: NewModel(cali, t, cf, 1),
		q:     &EMatrix{},
	}
	m.prop[0] = 1

	return
}

func (m *M0) SetDefault() {
	m.kappa = 1
	m.omega = 1
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
	*Model
	q0, q1, q2     *EMatrix
	kappa          float64
	omega0, omega2 float64
}

func NewH1(cali CodonSequences, t *tree.Tree, cf CodonFrequency) (m *H1) {
	m = &H1{
		Model: NewModel(cali, t, cf, 4),
		q0:    &EMatrix{},
		q1:    &EMatrix{},
		q2:    &EMatrix{},
	}
	return

}

func (m *H1) SetParameters(kappa float64, omega0, omega2 float64, p0, p1 float64) {
	m.kappa = kappa
	m.omega0 = omega0
	m.omega2 = omega2
	m.prop[0] = p0
	m.prop[1] = p1
	m.prop[2] = (1 - p0 - p1) * p0 / (p0 + p1)
	m.prop[3] = (1 - p0 - p1) * p1 / (p0 + p1)
	m.UpdateMatrices()
}

func (m *H1) SetDefault() {
	m.kappa = 1
	m.omega0 = 0.5
	m.omega2 = 2
	m.prop[0] = 0.25
	m.prop[1] = 0.25
	m.prop[2] = 0.25
	m.prop[3] = 0.25
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

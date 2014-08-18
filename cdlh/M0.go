package main

import (
	"math"

	"bitbucket.com/Davydov/golh/tree"
)

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

func (m *M0) GetNumberOfParameters() int {
	// root branch is not considered
	return 2 + m.tree.NNodes() - 1
}

func (m *M0) GetParameter(i int) float64 {
	switch i {
	case 0:
		return m.kappa
	case 1:
		return m.omega
	default:
		return m.tree.Nodes()[i-2+1].BranchLength
	}
}

func (m *M0) SetParameter(i int, value float64) {
	switch i {
	case 0:
		m.kappa = math.Abs(value)
		m.UpdateMatrix()
		m.ExpBranches()
	case 1:
		m.omega = math.Abs(value)
		m.UpdateMatrix()
		m.ExpBranches()
	default:
		br := i - 2 + 1
		m.tree.Nodes()[br].BranchLength = math.Abs(value)
		m.ExpBranch(br)
	}
}

func (m *M0) SetParameters(kappa, omega float64) {
	m.kappa = kappa
	m.omega = omega
	m.UpdateMatrix()
	m.ExpBranches()
}

func (m *M0) UpdateMatrix() {
	Q, s := createTransitionMatrix(m.cf, m.kappa, m.omega, m.q.Q)
	m.q.Set(Q, s)

	err := m.q.Eigen()
	if err != nil {
		panic("error finding eigen")
	}
	for i := 0; i < len(m.qs[0]); i++ {
		m.qs[0][i] = m.q
		m.scale[i] = m.q.Scale
	}
}

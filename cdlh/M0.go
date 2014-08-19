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

func (m *M0) SetDefaults() {
	m.kappa = 1
	m.omega = 1
	m.UpdateMatrix()
	m.ExpBranches()
}

func (m *M0) GetNumberOfParameters() (np int) {
	// root branch is not considered
	np = 2
	np += m.Model.GetNumberOfParameters()
	return
}

func (m *M0) GetParameterName(i int) string {
	switch i {
	case 0:
		return "kappa"
	case 1:
		return "omega"
	default:
		return m.Model.GetParameterName(i - 2)
	}
}

func (m *M0) GetParameter(i int) float64 {
	switch i {
	case 0:
		return m.kappa
	case 1:
		return m.omega
	default:
		return m.Model.GetParameter(i - 2)
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
		m.Model.SetParameter(i - 2, value)
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

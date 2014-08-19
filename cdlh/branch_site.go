package main

import (
	"math"

	"bitbucket.com/Davydov/golh/tree"
)

type BranchSite struct {
	*Model
	q0, q1, q2     *EMatrix
	kappa          float64
	omega0, omega2 float64
	p01sum, p0prop float64
}

func NewBranchSite(cali CodonSequences, t *tree.Tree, cf CodonFrequency) (m *BranchSite) {
	m = &BranchSite{
		Model: NewModel(cali, t, cf, 4),
		q0:    &EMatrix{},
		q1:    &EMatrix{},
		q2:    &EMatrix{},
	}
	return

}

func (m *BranchSite) SetParameters(kappa float64, omega0, omega2 float64, p0, p1 float64) {
	m.kappa = kappa
	m.omega0 = omega0
	m.omega2 = omega2
	m.p01sum = p0 + p1
	m.p0prop = p0 / (p0 + p1)
	m.UpdateMatrices(true, true, true)
	m.SetBranchMatrices()
	m.ExpBranches()
}

func (m *BranchSite) SetDefaults() {
	m.kappa = 1
	m.omega0 = 0.5
	m.omega2 = 2
	m.p01sum = 0.5
	m.p0prop = 0.5
	m.SetBranchMatrices()
	m.UpdateMatrices(true, true, true)
	m.ExpBranches()
}

func (m *BranchSite) GetNumberOfParameters() (np int) {
	// root branch is not considered
	np = 5
	np += m.Model.GetNumberOfParameters()
	return
}

func (m *BranchSite) GetParameterName(i int) string {
	switch i {
	case 0:
		return "kappa"
	case 1:
		return "omega0"
	case 2:
		return "omega2"
	case 3:
		// We use reparametrization
		return "p01sum"
	case 4:
		return "p0prop"
	default:
		return m.Model.GetParameterName(i - 5)
	}
}

func (m *BranchSite) GetParameter(i int) float64 {
	switch i {
	case 0:
		return m.kappa
	case 1:
		return m.omega0
	case 2:
		return m.omega2
	case 3:
		// We use reparametrization
		return m.p01sum
	case 4:
		return m.p0prop
	default:
		return m.Model.GetParameter(i - 5)
	}
}

func (m *BranchSite) SetParameter(i int, value float64) {
	switch i {
	case 0:
		m.kappa = math.Abs(value)
		m.UpdateMatrices(true, true, true)
		m.ExpBranches()
	case 1:
		m.omega0 = math.Abs(value)
		m.UpdateMatrices(true, false, false)
		m.ExpBranches()
	case 2:
		m.omega2 = Reflect(value, 1, math.Inf(+1))
		m.UpdateMatrices(false, false, true)
		m.ExpBranches()
	case 3:
		m.p01sum = Reflect(value, 1e-6, 1)
		m.UpdateMatrices(false, false, false)
	case 4:
		m.p0prop = Reflect(value, 0, 1)
		m.UpdateMatrices(false, false, false)
	default:
		m.Model.SetParameter(i - 5, value)
	}
}

func (m *BranchSite) SetBranchMatrices() {
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

func (m *BranchSite) UpdateMatrices(uq0, uq1, uq2 bool) {
	if uq0 {
		Q0, s0 := createTransitionMatrix(m.cf, m.kappa, m.omega0, m.q0.Q)
		m.q0.Set(Q0, s0)
		err := m.q0.Eigen()
		if err != nil {
			panic("error eigen q0")
		}
	}

	if uq1 {
		Q1, s1 := createTransitionMatrix(m.cf, m.kappa, 1, m.q1.Q)
		m.q1.Set(Q1, s1)
		err := m.q1.Eigen()
		if err != nil {
			panic("error eigen q1")
		}
	}

	if uq2 {
		Q2, s2 := createTransitionMatrix(m.cf, m.kappa, m.omega2, m.q2.Q)
		m.q2.Set(Q2, s2)
		err := m.q2.Eigen()
		if err != nil {
			panic("error eigen q2")
		}
	}

	p0 := m.p0prop * m.p01sum
	p1 := m.p01sum - p0
	m.prop[0] = p0
	m.prop[1] = p1
	m.prop[2] = (1 - p0 - p1) * p0 / (p0 + p1)
	m.prop[3] = (1 - p0 - p1) * p1 / (p0 + p1)

	for _, node := range m.tree.Nodes() {
		if node.Class == 0 {
			m.scale[node.Id] = (m.prop[0]+m.prop[2])*m.q0.Scale + (m.prop[1]+m.prop[3])*m.q1.Scale
		} else {
			m.scale[node.Id] = m.prop[0]*m.q0.Scale + m.prop[1]*m.q1.Scale + (m.prop[2]+m.prop[3])*m.q2.Scale
		}
	}

}

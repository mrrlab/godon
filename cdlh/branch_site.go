package main

import (
	"bitbucket.com/Davydov/golh/tree"
)

type BranchSite struct {
	*Model
	q0, q1, q2     *EMatrix
	kappa          float64
	omega0, omega2 float64
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
	m.prop[0] = p0
	m.prop[1] = p1
	m.prop[2] = (1 - p0 - p1) * p0 / (p0 + p1)
	m.prop[3] = (1 - p0 - p1) * p1 / (p0 + p1)
	m.UpdateMatrices(true, true, true)
	m.ExpBranches()
}

func (m *BranchSite) SetDefault() {
	m.kappa = 1
	m.omega0 = 0.5
	m.omega2 = 2
	m.prop[0] = 0.25
	m.prop[1] = 0.25
	m.prop[2] = 0.25
	m.prop[3] = 0.25
	m.UpdateMatrices(true, true, true)
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

	for _, node := range m.tree.Nodes() {
		if node.Class == 0 {
			m.scale[node.Id] = (m.prop[0]+m.prop[2])*m.q0.Scale + (m.prop[1]+m.prop[3])*m.q1.Scale
		} else {
			m.scale[node.Id] = m.prop[0]*m.q0.Scale + m.prop[1]*m.q1.Scale + (m.prop[2]+m.prop[3])*m.q2.Scale
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

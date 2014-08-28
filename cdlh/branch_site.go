package main

import (
	"bitbucket.com/Davydov/golh/mcmc"
	"bitbucket.com/Davydov/golh/tree"
)

type BranchSite struct {
	*Model
	q0, q1, q2     *EMatrix
	kappa          float64
	omega0, omega2 float64
	p01sum, p0prop float64
	parameters     mcmc.Parameters
	q0done, q1done, q2done bool
	propdone bool
}

func NewBranchSite(cali CodonSequences, t *tree.Tree, cf CodonFrequency, optBranch bool) (m *BranchSite) {
	m = &BranchSite{
		Model: NewModel(cali, t, cf, 4, optBranch),
		q0:    &EMatrix{},
		q1:    &EMatrix{},
		q2:    &EMatrix{},
	}
	m.parameters = m.Model.parameters
	kappa := mcmc.NewFloat64Parameter(&m.kappa, "kappa")
	kappa.OnChange = func() {
		m.q0done = false
		m.q1done = false
		m.q2done = false
		m.expAllBr = false
	}
	kappa.PriorFunc = mcmc.UniformPrior(0, 20, false, true)
	kappa.ProposalFunc = mcmc.NormalProposal(0.01)
	kappa.Min = 0
	kappa.Max = 20
	m.parameters = append(m.parameters, kappa)

	omega0 := mcmc.NewFloat64Parameter(&m.omega0, "omega0")
	omega0.OnChange = func() {
		m.q0done = false
		m.expAllBr = false
	}
	omega0.PriorFunc = mcmc.GammaPrior(1, 2, false)
	omega0.ProposalFunc = mcmc.NormalProposal(0.01)
	omega0.Min = 0
	m.parameters = append(m.parameters, omega0)

	omega2 := mcmc.NewFloat64Parameter(&m.omega2, "omega2")
	omega2.OnChange = func() {
		m.q2done = false
		m.expAllBr = false
	}
	omega2.PriorFunc = mcmc.GammaPrior(1, 2, false)
	omega2.ProposalFunc = mcmc.NormalProposal(0.01)
	omega2.Min = 0
	m.parameters = append(m.parameters, omega2)

	p01sum := mcmc.NewFloat64Parameter(&m.p01sum, "p01sum")
	p01sum.OnChange = func() {
		m.propdone = false
	}
	p01sum.PriorFunc = mcmc.UniformPrior(0, 1, false, false)
	p01sum.ProposalFunc = mcmc.NormalProposal(0.01)
	m.parameters = append(m.parameters, p01sum)

	p0prop := mcmc.NewFloat64Parameter(&m.p0prop, "p0prop")
	p0prop.OnChange = func() {
		m.propdone = false
	}
	p0prop.PriorFunc = mcmc.UniformPrior(0, 1, false, false)
	p0prop.ProposalFunc = mcmc.NormalProposal(0.01)
	m.parameters = append(m.parameters, p0prop)

	m.SetBranchMatrices()
	m.SetDefaults()

	return

}

func (m *BranchSite) GetParameters() mcmc.Parameters {
	return m.parameters
}

func (m *BranchSite) SetParameters(kappa float64, omega0, omega2 float64, p0, p1 float64) {
	m.kappa = kappa
	m.omega0 = omega0
	m.omega2 = omega2
	m.p01sum = p0 + p1
	m.p0prop = p0 / (p0 + p1)
}

func (m *BranchSite) SetDefaults() {
	m.SetParameters(1, 0.5, 2, 0.5, 0.5)
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

func (m *BranchSite) UpdateProportions() {
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
	m.propdone = true
}

func (m *BranchSite) UpdateMatrices() {
	if !m.q0done {
		Q0, s0 := createTransitionMatrix(m.cf, m.kappa, m.omega0, m.q0.Q)
		m.q0.Set(Q0, s0)
		err := m.q0.Eigen()
		if err != nil {
			panic("error eigen q0")
		}
		m.q0done = true
	}

	if !m.q1done {
		Q1, s1 := createTransitionMatrix(m.cf, m.kappa, 1, m.q1.Q)
		m.q1.Set(Q1, s1)
		err := m.q1.Eigen()
		if err != nil {
			panic("error eigen q1")
		}
		m.q1done = true
	}

	if !m.q2done {
		Q2, s2 := createTransitionMatrix(m.cf, m.kappa, m.omega2, m.q2.Q)
		m.q2.Set(Q2, s2)
		err := m.q2.Eigen()
		if err != nil {
			panic("error eigen q2")
		}
		m.q2done = true
	}

		m.UpdateProportions()
}

func (m *BranchSite) Likelihood() float64 {
	if !m.q0done || !m.q1done || !m.q2done {
		m.UpdateMatrices()
	}
	if !m.propdone {
		m.UpdateProportions()
	}
	return m.Model.Likelihood()
}

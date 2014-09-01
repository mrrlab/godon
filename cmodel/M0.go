package cmodel

import (
	"bitbucket.com/Davydov/golh/mcmc"
	"bitbucket.com/Davydov/golh/tree"
)

type M0 struct {
	*Model
	q            *EMatrix
	omega, kappa float64
	parameters   mcmc.Parameters
	qdone        bool
}

func NewM0(cali CodonSequences, t *tree.Tree, cf CodonFrequency, optBranch bool) (m *M0) {
	m = &M0{
		Model: NewModel(cali, t, cf, 1, optBranch),
		q:     &EMatrix{},
	}
	m.prop[0] = 1
	m.parameters = m.Model.parameters

	m.addParameters()
	m.SetDefaults()
	return
}

func (m *M0) SetAdaptive(as *mcmc.AdaptiveSettings) {
	m.Model.SetAdaptive()
	m.parameters = m.Model.parameters
	m.addAdaptiveParameters(as)
}

func (m *M0) addParameters() {
	omega := mcmc.NewFloat64Parameter(&m.omega, "omega")
	omega.OnChange = func() {
		m.qdone = false
		m.expAllBr = false
	}
	omega.PriorFunc = mcmc.GammaPrior(1, 2, false)
	omega.ProposalFunc = mcmc.NormalProposal(0.01)
	omega.Min = 0

	kappa := mcmc.NewFloat64Parameter(&m.kappa, "kappa")
	kappa.OnChange = func() {
		m.qdone = false
		m.expAllBr = false
	}
	kappa.PriorFunc = mcmc.UniformPrior(0, 20, false, true)
	kappa.ProposalFunc = mcmc.NormalProposal(0.01)
	kappa.Min = 0
	kappa.Max = 20

	m.parameters = append(m.parameters, omega)
	m.parameters = append(m.parameters, kappa)
}

func (m *M0) addAdaptiveParameters(as *mcmc.AdaptiveSettings) {
	omega := mcmc.NewAdaptiveParameter(&m.omega, "omega", as)
	omega.OnChange = func() {
		m.qdone = false
		m.expAllBr = false
	}
	omega.PriorFunc = mcmc.GammaPrior(1, 2, false)
	omega.Min = 0

	kappa := mcmc.NewAdaptiveParameter(&m.kappa, "kappa", as)
	kappa.OnChange = func() {
		m.qdone = false
		m.expAllBr = false
	}
	kappa.PriorFunc = mcmc.UniformPrior(0, 20, false, true)
	kappa.ProposalFunc = mcmc.NormalProposal(0.01)
	kappa.Min = 0
	kappa.Max = 20

	m.parameters = append(m.parameters, omega)
	m.parameters = append(m.parameters, kappa)
}

func (m *M0) GetModelParameters() mcmc.Parameters {
	return mcmc.Parameters(m.parameters)
}

func (m *M0) GetParameters() (kappa, omega float64) {
	return m.kappa, m.omega
}

func (m *M0) SetParameters(kappa, omega float64) {
	m.kappa = kappa
	m.omega = omega
	m.qdone = false
}

func (m *M0) SetDefaults() {
	m.SetParameters(1, 1)
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
	m.expAllBr = false
}

func (m *M0) Likelihood() float64 {
	if !m.qdone {
		m.UpdateMatrix()
	}
	return m.Model.Likelihood()
}

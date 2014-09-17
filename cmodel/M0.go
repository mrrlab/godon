package cmodel

import (
	"bitbucket.com/Davydov/golh/optimize"
	"bitbucket.com/Davydov/golh/tree"
)

type M0 struct {
	*Model
	q            *EMatrix
	omega, kappa float64
	parameters   optimize.FloatParameters
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

func (m *M0) Copy() optimize.Optimizable {
	newM := &M0{
		Model: NewModel(m.cali, m.tree.Copy(), m.cf, 1, m.optBranch),
		q:     &EMatrix{},
		omega: m.omega,
		kappa: m.kappa,
	}
	newM.prop[0] = 1
	newM.as = m.as
	newM.Model.setParameters()
	newM.parameters = newM.Model.parameters

	if m.as != nil {
		newM.addAdaptiveParameters()
	} else {
		newM.addParameters()
	}
	return newM
}

func (m *M0) SetAdaptive(as *optimize.AdaptiveSettings) {
	m.as = as
	m.Model.setParameters()
	m.parameters = m.Model.parameters
	m.addAdaptiveParameters()
}

func (m *M0) addParameters() {
	omega := optimize.NewBasicFloatParameter(&m.omega, "omega")
	omega.OnChange = func() {
		m.qdone = false
		m.expAllBr = false
	}
	omega.PriorFunc = optimize.GammaPrior(1, 2, false)
	omega.ProposalFunc = optimize.NormalProposal(0.01)
	omega.Min = 0

	kappa := optimize.NewBasicFloatParameter(&m.kappa, "kappa")
	kappa.OnChange = func() {
		m.qdone = false
		m.expAllBr = false
	}
	kappa.PriorFunc = optimize.UniformPrior(0, 20, false, true)
	kappa.ProposalFunc = optimize.NormalProposal(0.01)
	kappa.Min = 0
	kappa.Max = 20

	m.parameters.Append(omega)
	m.parameters.Append(kappa)
}

func (m *M0) addAdaptiveParameters() {
	omega := optimize.NewAdaptiveParameter(&m.omega, "omega", m.as)
	omega.OnChange = func() {
		m.qdone = false
		m.expAllBr = false
	}
	omega.PriorFunc = optimize.GammaPrior(1, 2, false)
	omega.Min = 0

	kappa := optimize.NewAdaptiveParameter(&m.kappa, "kappa", m.as)
	kappa.OnChange = func() {
		m.qdone = false
		m.expAllBr = false
	}
	kappa.PriorFunc = optimize.UniformPrior(0, 20, false, true)
	kappa.ProposalFunc = optimize.NormalProposal(0.01)
	kappa.Min = 0
	kappa.Max = 20

	m.parameters.Append(omega)
	m.parameters.Append(kappa)
}

func (m *M0) GetFloatParameters() optimize.FloatParameters {
	return m.parameters
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

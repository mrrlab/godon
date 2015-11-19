package cmodel

import (
	"bitbucket.com/Davydov/godon/optimize"
	"bitbucket.com/Davydov/godon/tree"
)

type M0vrate struct {
	*Model
	q0, q1       *EMatrix
	omega, kappa float64
	p            float64 // proportion of non-scaled
	s            float64 // scale-factor
	parameters   optimize.FloatParameters
	qdone        bool
	propdone     bool
}

func NewM0vrate(cali CodonSequences, t *tree.Tree, cf CodonFrequency, optBranch bool) (m *M0vrate) {
	m = &M0vrate{
		Model: NewModel(cali, t, cf, 2, optBranch),
		q0:    &EMatrix{},
		q1:    &EMatrix{},
	}
	m.parameters = m.Model.parameters

	m.addParameters()
	m.SetDefaults()
	return
}

func (m *M0vrate) Copy() optimize.Optimizable {
	newM := &M0vrate{
		Model: NewModel(m.cali, m.tree.Copy(), m.cf, 2, m.optBranch),
		q0:    &EMatrix{},
		q1:    &EMatrix{},
		omega: m.omega,
		kappa: m.kappa,
		p:     m.p,
		s:     m.s,
	}
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

func (m *M0vrate) SetAdaptive(as *optimize.AdaptiveSettings) {
	m.as = as
	m.Model.setParameters()
	m.parameters = m.Model.parameters
	m.addAdaptiveParameters()
}

func (m *M0vrate) addParameters() {
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

	p := optimize.NewBasicFloatParameter(&m.p, "p")
	p.OnChange = func() {
		m.propdone = false
	}
	p.PriorFunc = optimize.UniformPrior(0, 1, true, true)
	p.ProposalFunc = optimize.NormalProposal(0.01)
	p.Min = 0
	p.Max = 1

	s := optimize.NewBasicFloatParameter(&m.s, "s")
	s.OnChange = func() {
		m.qdone = false
		m.expAllBr = false
	}
	s.PriorFunc = optimize.ExponentialPrior(1, false)
	s.ProposalFunc = optimize.NormalProposal(0.01)
	s.Min = 1. / 100
	s.Max = 100

	m.parameters.Append(omega)
	m.parameters.Append(kappa)
	m.parameters.Append(p)
	m.parameters.Append(s)
}

func (m *M0vrate) addAdaptiveParameters() {
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

	p := optimize.NewAdaptiveParameter(&m.p, "p", m.as)
	p.OnChange = func() {
		m.propdone = false
	}
	p.PriorFunc = optimize.UniformPrior(0, 1, true, true)
	p.ProposalFunc = optimize.NormalProposal(0.01)
	p.Min = 0
	p.Max = 1

	s := optimize.NewAdaptiveParameter(&m.s, "s", m.as)
	s.OnChange = func() {
		m.qdone = false
		m.expAllBr = false
	}
	s.PriorFunc = optimize.ExponentialPrior(1, false)
	s.ProposalFunc = optimize.NormalProposal(0.01)
	s.Min = 1. / 100
	s.Max = 100

	m.parameters.Append(omega)
	m.parameters.Append(kappa)
	m.parameters.Append(p)
	m.parameters.Append(s)
}

func (m *M0vrate) GetFloatParameters() optimize.FloatParameters {
	return m.parameters
}

func (m *M0vrate) GetParameters() (kappa, omega float64) {
	return m.kappa, m.omega
}

func (m *M0vrate) SetParameters(kappa, omega, p, s float64) {
	m.kappa = kappa
	m.omega = omega
	m.p = p
	m.s = s
	m.qdone = false
}

func (m *M0vrate) SetDefaults() {
	m.SetParameters(1, 1, 0.5, 2)
}

func (m *M0vrate) UpdateProportions() {
	m.prop[0] = m.p
	m.prop[1] = 1 - m.p

	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		m.scale[node.Id] = m.prop[0]*m.q0.Scale + m.prop[1]*m.q1.Scale
		m.qs[0][node.Id] = m.q0
		m.qs[1][node.Id] = m.q1
	}
	m.expAllBr = false
}

func (m *M0vrate) UpdateMatrix() {
	Q0, scale := createTransitionMatrix(m.cf, m.kappa, m.omega, m.q0.Q)
	m.q0.Set(Q0, scale)
	err := m.q0.Eigen()
	if err != nil {
		panic("error finding eigen")
	}

	Q1, scale := createTransitionMatrix(m.cf, m.kappa, m.omega, m.q1.Q)
	m.q1.Set(Q1, scale/m.s)
	err = m.q1.Eigen()
	if err != nil {
		panic("error finding eigen")
	}

	m.UpdateProportions()
}

func (m *M0vrate) Likelihood() float64 {
	if !m.qdone {
		m.UpdateMatrix()
	}
	m.UpdateProportions()
	return m.Model.Likelihood()
}

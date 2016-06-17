package cmodel

import (
	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/tree"
)

type M0vrate struct {
	*BaseModel
	q0, q1       *codon.EMatrix
	omega, kappa float64
	p            float64 // proportion of non-scaled
	s            float64 // scale-factor
	qdone        bool
	propdone     bool
}

func NewM0vrate(cali codon.CodonSequences, t *tree.Tree, cf codon.CodonFrequency) (m *M0vrate) {
	m = &M0vrate{
		q0: &codon.EMatrix{CF: cf},
		q1: &codon.EMatrix{CF: cf},
	}

	m.BaseModel = NewBaseModel(cali, t, cf, m)

	m.setupParameters()
	m.SetDefaults()
	return
}

func (m *M0vrate) GetNClass() int {
	return 2
}

func (m *M0vrate) Copy() optimize.Optimizable {
	newM := &M0vrate{
		BaseModel: m.BaseModel.Copy(),
		q0:        &codon.EMatrix{CF: m.cf},
		q1:        &codon.EMatrix{CF: m.cf},
		omega:     m.omega,
		kappa:     m.kappa,
		s:         m.s,
	}
	newM.BaseModel.Model = newM

	newM.setupParameters()
	return newM
}

func (m *M0vrate) addParameters(fpg optimize.FloatParameterGenerator) {
	omega := fpg(&m.omega, "omega")
	omega.SetOnChange(func() {
		m.qdone = false
		m.expAllBr = false
	})
	omega.SetPriorFunc(optimize.GammaPrior(1, 2, false))
	omega.SetProposalFunc(optimize.NormalProposal(0.01))
	omega.SetMin(1)
	omega.SetMax(1000)

	kappa := fpg(&m.kappa, "kappa")
	kappa.SetOnChange(func() {
		m.qdone = false
		m.expAllBr = false
	})
	kappa.SetPriorFunc(optimize.UniformPrior(0, 20, false, true))
	kappa.SetProposalFunc(optimize.NormalProposal(0.01))
	kappa.SetMin(1e-2)
	kappa.SetMax(20)

	s := fpg(&m.s, "s")
	s.SetOnChange(func() {
		m.qdone = false
		m.expAllBr = false
	})
	s.SetPriorFunc(optimize.ExponentialPrior(1, false))
	s.SetProposalFunc(optimize.NormalProposal(0.01))
	s.SetMin(1. / 1000)
	s.SetMax(100)

	m.parameters.Append(omega)
	m.parameters.Append(kappa)
	m.parameters.Append(s)
}

func (m *M0vrate) GetParameters() (kappa, omega float64) {
	return m.kappa, m.omega
}

func (m *M0vrate) SetParameters(kappa, omega, s float64) {
	m.kappa = kappa
	m.omega = omega
	m.s = s
	m.qdone = false
}

func (m *M0vrate) SetDefaults() {
	m.SetParameters(1, 1, 2)
}

func (m *M0vrate) updateProportions() {
	m.prop[0] = 1 / m.s
	m.prop[1] = 1 - m.prop[0]

	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		m.scale[node.Id] = m.prop[0]*m.q0.Scale + m.prop[1]*m.q1.Scale
		m.qs[0][node.Id] = m.q0
		m.qs[1][node.Id] = m.q1
	}
	m.propdone = true
	m.expAllBr = false
}

func (m *M0vrate) updateMatrix() {
	Q0, scale := codon.CreateTransitionMatrix(m.cf, m.kappa, m.omega, m.q0.Q)
	m.q0.Set(Q0, scale)
	err := m.q0.Eigen()
	if err != nil {
		panic("error finding eigen")
	}
	m.q1 = m.q0.Copy(nil)
	m.q1.ScaleD(m.s)

	m.updateProportions()
}

func (m *M0vrate) Likelihood() float64 {
	if !m.qdone {
		m.updateMatrix()
	}
	if !m.propdone {
		m.updateProportions()
	}
	return m.BaseModel.Likelihood()
}

package cmodel

import (
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/tree"
)

type BranchSiteC struct {
	*BaseModel
	q0, q1, q2             *EMatrix
	kappa                  float64
	omega0, omega2         float64
	p0prop                 float64
	q0done, q1done, q2done bool
	propdone               bool
}

func NewBranchSiteC(cali CodonSequences, t *tree.Tree, cf CodonFrequency) (m *BranchSiteC) {
	m = &BranchSiteC{
		q0: &EMatrix{},
		q1: &EMatrix{},
		q2: &EMatrix{},
	}

	m.BaseModel = NewBaseModel(cali, t, cf, m)

	m.addParameters()
	m.SetBranchMatrices()
	m.SetDefaults()

	return
}

func (m *BranchSiteC) GetNClass() int {
	return 4
}

func (m *BranchSiteC) Copy() optimize.Optimizable {
	newM := &BranchSiteC{
		BaseModel: m.BaseModel.Copy(),
		q0:        &EMatrix{},
		q1:        &EMatrix{},
		q2:        &EMatrix{},
		kappa:     m.kappa,
		omega0:    m.omega0,
		omega2:    m.omega2,
		p0prop:    m.p0prop,
	}
	newM.BaseModel.Model = newM

	newM.setupParameters()
	newM.SetBranchMatrices()
	return newM
}

func (m *BranchSiteC) addParameters() {
	kappa := optimize.NewBasicFloatParameter(&m.kappa, "kappa")
	kappa.OnChange = func() {
		m.q0done = false
		m.q1done = false
		m.q2done = false
	}
	kappa.PriorFunc = optimize.UniformPrior(0, 20, false, true)
	kappa.ProposalFunc = optimize.NormalProposal(0.01)
	kappa.Min = 0
	kappa.Max = 20
	m.parameters = append(m.parameters, kappa)

	omega0 := optimize.NewBasicFloatParameter(&m.omega0, "omega0")
	omega0.OnChange = func() {
		m.q0done = false
	}
	omega0.PriorFunc = optimize.GammaPrior(1, 2, false)
	omega0.ProposalFunc = optimize.NormalProposal(0.01)
	omega0.Min = 0
	m.parameters = append(m.parameters, omega0)

	omega2 := optimize.NewBasicFloatParameter(&m.omega2, "omega2")
	omega2.OnChange = func() {
		m.q2done = false
	}
	omega2.PriorFunc = optimize.GammaPrior(1, 2, false)
	omega2.ProposalFunc = optimize.NormalProposal(0.01)
	omega2.Min = 1
	m.parameters = append(m.parameters, omega2)

	p0prop := optimize.NewBasicFloatParameter(&m.p0prop, "p0prop")
	p0prop.OnChange = func() {
		m.propdone = false
	}
	p0prop.PriorFunc = optimize.UniformPrior(0, 1, false, false)
	p0prop.Min = 0
	p0prop.Max = 1
	p0prop.ProposalFunc = optimize.NormalProposal(0.01)
	m.parameters = append(m.parameters, p0prop)
}

func (m *BranchSiteC) addAdaptiveParameters() {
	kappa := optimize.NewAdaptiveParameter(&m.kappa, "kappa", m.as)
	kappa.OnChange = func() {
		m.q0done = false
		m.q1done = false
		m.q2done = false
	}
	kappa.PriorFunc = optimize.UniformPrior(0, 20, false, true)
	kappa.Min = 0
	kappa.Max = 20
	m.parameters = append(m.parameters, kappa)

	omega0 := optimize.NewAdaptiveParameter(&m.omega0, "omega0", m.as)
	omega0.OnChange = func() {
		m.q0done = false
	}
	omega0.PriorFunc = optimize.GammaPrior(1, 2, false)
	omega0.Min = 0
	m.parameters = append(m.parameters, omega0)

	omega2 := optimize.NewAdaptiveParameter(&m.omega2, "omega2", m.as)
	omega2.OnChange = func() {
		m.q2done = false
	}
	omega2.PriorFunc = optimize.GammaPrior(1, 2, false)
	omega2.Min = 1
	m.parameters = append(m.parameters, omega2)

	p0prop := optimize.NewAdaptiveParameter(&m.p0prop, "p0prop", m.as)
	p0prop.Min = 0
	p0prop.Max = 1
	p0prop.OnChange = func() {
		m.propdone = false
	}
	p0prop.PriorFunc = optimize.UniformPrior(0, 1, false, false)
	m.parameters = append(m.parameters, p0prop)
}

func (m *BranchSiteC) SetParameters(kappa float64, omega0, omega2 float64, p0prop float64) {
	m.kappa = kappa
	m.omega0 = omega0

	m.omega2 = omega2

	m.p0prop = p0prop

	m.q0done, m.q1done, m.q2done = false, false, false
}

func (m *BranchSiteC) GetParameters() (kappa float64, omega0, omega2 float64, p0 float64) {
	return m.kappa, m.omega0, m.omega2, m.p0prop
}

func (m *BranchSiteC) SetDefaults() {
	m.SetParameters(1, 0.5, 2, 0.5)
}

func (m *BranchSiteC) SetBranchMatrices() {
	for i := 0; i < len(m.qs); i++ {
		for _, node := range m.tree.NodeIdArray() {
			if node == nil {
				continue
			}
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

func (m *BranchSiteC) UpdateProportions() {
	p2 := m.omega2 / 1000
	m.prop[0] = m.p0prop * (1 - p2)
	m.prop[1] = (1 - m.p0prop) * (1 - p2)
	m.prop[2] = m.p0prop * p2
	m.prop[3] = (1 - m.p0prop) * p2

	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		if node.Class == 0 {
			m.scale[node.Id] = (m.prop[0]+m.prop[2])*m.q0.Scale + (m.prop[1]+m.prop[3])*m.q1.Scale
		} else {
			m.scale[node.Id] = m.prop[0]*m.q0.Scale + m.prop[1]*m.q1.Scale + (m.prop[2]+m.prop[3])*m.q2.Scale
		}
	}
	m.propdone = true
	m.expAllBr = false
}

func (m *BranchSiteC) UpdateMatrices() {
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
	m.expAllBr = false
}

func (m *BranchSiteC) Likelihood() float64 {
	if !m.q0done || !m.q1done || !m.q2done {
		m.UpdateMatrices()
	}
	if !m.propdone {
		m.UpdateProportions()
	}
	return m.BaseModel.Likelihood()
}

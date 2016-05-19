package cmodel

import (
	"math"
	"math/rand"

	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/tree"
)

type BranchSite struct {
	*BaseModel
	q0, q1, q2             *EMatrix
	kappa                  float64
	omega0, omega2         float64
	p01sum, p0prop         float64
	fixw2                  bool
	q0done, q1done, q2done bool
	propdone               bool
}

func NewBranchSite(cali CodonSequences, t *tree.Tree, cf CodonFrequency, fixw2 bool) (m *BranchSite) {
	m = &BranchSite{
		fixw2: fixw2,
		q0:    &EMatrix{},
		q1:    &EMatrix{},
		q2:    &EMatrix{},
	}

	m.BaseModel = NewBaseModel(cali, t, cf, m)

	m.setupParameters()
	m.setBranchMatrices()
	m.SetDefaults()

	return

}

func (m *BranchSite) GetNClass() int {
	return 4
}

func (m *BranchSite) Copy() optimize.Optimizable {
	newM := &BranchSite{
		BaseModel: m.BaseModel.Copy(),
		q0:        &EMatrix{},
		q1:        &EMatrix{},
		q2:        &EMatrix{},
		kappa:     m.kappa,
		omega0:    m.omega0,
		omega2:    m.omega2,
		p01sum:    m.p01sum,
		p0prop:    m.p0prop,
		fixw2:     m.fixw2,
	}
	newM.BaseModel.Model = newM
	newM.setupParameters()
	newM.setBranchMatrices()
	return newM
}

func (m *BranchSite) addParameters(fpg optimize.FloatParameterGenerator) {
	kappa := fpg(&m.kappa, "kappa")
	kappa.SetOnChange(func() {
		m.q0done = false
		m.q1done = false
		m.q2done = false
	})
	kappa.SetPriorFunc(optimize.UniformPrior(0, 20, false, true))
	kappa.SetProposalFunc(optimize.NormalProposal(0.01))
	kappa.SetMin(0.1)
	kappa.SetMax(20)
	m.parameters.Append(kappa)

	omega0 := fpg(&m.omega0, "omega0")
	omega0.SetOnChange(func() {
		m.q0done = false
	})
	omega0.SetPriorFunc(optimize.GammaPrior(1, 2, false))
	omega0.SetProposalFunc(optimize.NormalProposal(0.01))
	omega0.SetMin(0)
	omega0.SetMax(1)
	m.parameters.Append(omega0)

	if !m.fixw2 {
		omega2 := fpg(&m.omega2, "omega2")
		omega2.SetOnChange(func() {
			m.q2done = false
		})
		omega2.SetPriorFunc(optimize.GammaPrior(1, 2, false))
		omega2.SetProposalFunc(optimize.NormalProposal(0.01))
		omega2.SetMin(1)
		m.parameters.Append(omega2)
	}

	p01sum := fpg(&m.p01sum, "p01sum")
	p01sum.SetOnChange(func() {
		m.propdone = false
	})
	p01sum.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
	p01sum.SetMin(0)
	p01sum.SetMax(1)
	p01sum.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(p01sum)

	p0prop := fpg(&m.p0prop, "p0prop")
	p0prop.SetOnChange(func() {
		m.propdone = false
	})
	p0prop.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
	p0prop.SetMin(0)
	p0prop.SetMax(1)
	p0prop.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(p0prop)
}

func (m *BranchSite) SetParameters(kappa float64, omega0, omega2 float64, p0, p1 float64) {
	m.kappa = kappa
	m.omega0 = omega0
	if !m.fixw2 {
		m.omega2 = omega2
	} else {
		m.omega2 = 1
	}
	m.p01sum = p0 + p1
	m.p0prop = p0 / (p0 + p1)
	m.q0done, m.q1done, m.q2done = false, false, false
}

func (m *BranchSite) GetParameters() (kappa float64, omega0, omega2 float64, p0, p1 float64) {
	return m.kappa, m.omega0, m.omega2, m.p01sum * m.p0prop, m.p01sum * (1 - m.p0prop)
}

func (m *BranchSite) SetDefaults() {
	// these parameters come from codeml
	kappa := 1.0
	omega0 := 0.2 + 0.1*rand.Float64()
	omega2 := 3.1 + rand.Float64()
	x0 := 1.0 + 0.5*rand.Float64()
	x1 := 0.2 * rand.Float64()
	p0 := math.Exp(x0) / (1 + math.Exp(x0) + math.Exp(x1))
	p1 := math.Exp(x1) / (1 + math.Exp(x0) + math.Exp(x1))
	m.SetParameters(kappa, omega0, omega2, p0, p1)
}

func (m *BranchSite) setBranchMatrices() {
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

func (m *BranchSite) updateProportions() {
	p0 := m.p0prop * m.p01sum
	p1 := m.p01sum - p0
	m.prop[0] = p0
	m.prop[1] = p1
	m.prop[2] = (1 - p0 - p1) * p0 / (p0 + p1)
	m.prop[3] = (1 - p0 - p1) * p1 / (p0 + p1)

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

func (m *BranchSite) updateMatrices() {
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

	m.updateProportions()
	m.expAllBr = false
}

func (m *BranchSite) Likelihood() float64 {
	if !m.q0done || !m.q1done || !m.q2done {
		m.updateMatrices()
	}
	if !m.propdone {
		m.updateProportions()
	}
	return m.BaseModel.Likelihood()
}

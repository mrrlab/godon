package cmodel

import (
	"math"
	"math/rand"

	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/paml"
	"bitbucket.org/Davydov/godon/tree"
)

type BranchSiteGamma struct {
	*BaseModel
	q0, q1, q2             *EMatrix
	q0s                    []*EMatrix
	q1s                    []*EMatrix
	q2s                    []*EMatrix
	kappa                  float64
	omega0, omega2         float64
	p01sum, p0prop         float64
	alpha                  float64
	rates                  []float64
	ncat                   int
	fixw2                  bool
	q0done, q1done, q2done bool
	propdone               bool
	ratesdone              bool
}

func NewBranchSiteGamma(cali CodonSequences, t *tree.Tree, cf CodonFrequency, ncat int, fixw2 bool) (m *BranchSiteGamma) {
	m = &BranchSiteGamma{
		fixw2: fixw2,
		ncat:  ncat,
		q0:    &EMatrix{},
		q1:    &EMatrix{},
		q2:    &EMatrix{},
		q0s:   make([]*EMatrix, ncat),
		q1s:   make([]*EMatrix, ncat),
		q2s:   make([]*EMatrix, ncat),
	}
	m.BaseModel = NewBaseModel(cali, t, cf, m)

	for i := 0; i < ncat; i++ {
		m.q0s[i] = &EMatrix{}
		m.q1s[i] = &EMatrix{}
		m.q2s[i] = &EMatrix{}
	}

	m.setupParameters()
	m.setBranchMatrices()
	m.SetDefaults()

	return

}

func (m *BranchSiteGamma) GetNClass() int {
	return 4 * m.ncat
}

func (m *BranchSiteGamma) Copy() optimize.Optimizable {
	newM := &BranchSiteGamma{
		BaseModel: m.BaseModel.Copy(),
		q0:        &EMatrix{},
		q1:        &EMatrix{},
		q2:        &EMatrix{},
		q0s:       make([]*EMatrix, m.ncat),
		q1s:       make([]*EMatrix, m.ncat),
		q2s:       make([]*EMatrix, m.ncat),
		kappa:     m.kappa,
		omega0:    m.omega0,
		omega2:    m.omega2,
		p01sum:    m.p01sum,
		p0prop:    m.p0prop,
		alpha:     m.alpha,
		ncat:      m.ncat,
		fixw2:     m.fixw2,
	}
	newM.BaseModel.Model = newM

	for i := 0; i < m.ncat; i++ {
		newM.q0s[i] = &EMatrix{}
		newM.q1s[i] = &EMatrix{}
		newM.q2s[i] = &EMatrix{}
	}

	newM.setupParameters()
	newM.setBranchMatrices()
	return newM
}

func (m *BranchSiteGamma) addParameters(fpg optimize.FloatParameterGenerator) {
	kappa := fpg(&m.kappa, "kappa")
	kappa.SetOnChange(func() {
		m.q0done = false
		m.q1done = false
		m.q2done = false
	})
	kappa.SetPriorFunc(optimize.UniformPrior(0, 20, false, true))
	kappa.SetProposalFunc(optimize.NormalProposal(0.01))
	kappa.SetMin(0)
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

	alpha := fpg(&m.alpha, "alpha")
	alpha.SetOnChange(func() {
		m.ratesdone = false
	})
	alpha.SetPriorFunc(optimize.GammaPrior(1, 2, false))
	alpha.SetMin(0)
	alpha.SetMax(1000)
	alpha.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(alpha)
}

func (m *BranchSiteGamma) SetParameters(kappa float64, omega0, omega2 float64, p0, p1 float64, alpha float64) {
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
	m.alpha = alpha
}

func (m *BranchSiteGamma) GetParameters() (kappa float64, omega0, omega2 float64, p0, p1 float64, alpha float64) {
	return m.kappa, m.omega0, m.omega2, m.p01sum * m.p0prop, m.p01sum * (1 - m.p0prop), m.alpha
}

func (m *BranchSiteGamma) SetDefaults() {
	// these parameters come from codeml
	kappa := 1.0
	omega0 := 0.2 + 0.1*rand.Float64()
	omega2 := 3.1 + rand.Float64()
	x0 := 1.0 + 0.5*rand.Float64()
	x1 := 0.2 * rand.Float64()
	p0 := math.Exp(x0) / (1 + math.Exp(x0) + math.Exp(x1))
	p1 := math.Exp(x1) / (1 + math.Exp(x0) + math.Exp(x1))
	alpha := 1.0
	m.SetParameters(kappa, omega0, omega2, p0, p1, alpha)
}

func (m *BranchSiteGamma) setBranchMatrices() {
	for i := 0; i < m.ncat; i++ {
		for j := 0; j < 4; j++ {
			for _, node := range m.tree.NodeIdArray() {
				if node == nil {
					continue
				}
				switch j {
				case 0:
					m.qs[i*4+j][node.Id] = m.q0s[i]
				case 1:
					m.qs[i*4+j][node.Id] = m.q1s[i]
				case 2:
					if node.Class == 0 {
						m.qs[i*4+j][node.Id] = m.q0s[i]
					} else {
						m.qs[i*4+j][node.Id] = m.q2s[i]
					}
				case 3:
					if node.Class == 0 {
						m.qs[i*4+j][node.Id] = m.q1s[i]
					} else {
						m.qs[i*4+j][node.Id] = m.q2s[i]
					}
				}
			}
		}
	}
}

func (m *BranchSiteGamma) updateProportions() {
	p0 := m.p0prop * m.p01sum
	p1 := m.p01sum - p0
	for i := 0; i < m.ncat; i++ {
		m.prop[i*4+0] = p0 / float64(m.ncat)
		m.prop[i*4+1] = p1 / float64(m.ncat)
		m.prop[i*4+2] = (1 - p0 - p1) * p0 / (p0 + p1) / float64(m.ncat)
		m.prop[i*4+3] = (1 - p0 - p1) * p1 / (p0 + p1) / float64(m.ncat)
	}

	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		if node.Class == 0 {
			m.scale[node.Id] = (m.prop[0]+m.prop[2])*float64(m.ncat)*m.q0.Scale + (m.prop[1]+m.prop[3])*float64(m.ncat)*m.q1.Scale
		} else {
			m.scale[node.Id] = m.prop[0]*float64(m.ncat)*m.q0.Scale + m.prop[1]*float64(m.ncat)*m.q1.Scale + (m.prop[2]+m.prop[3])*float64(m.ncat)*m.q2.Scale
		}
	}
	m.propdone = true
	m.expAllBr = false
}

func (m *BranchSiteGamma) updateMatrices() {
	if !m.q0done {
		Q0, s0 := createTransitionMatrix(m.cf, m.kappa, m.omega0, m.q0.Q)
		m.q0.Set(Q0, s0)
		err := m.q0.Eigen()
		if err != nil {
			panic("error eigen q0")
		}
		for i := 0; i < m.ncat; i++ {
			m.q0.Copy(m.q0s[i])
			m.q0s[i].ScaleD(m.rates[i])
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
		for i := 0; i < m.ncat; i++ {
			m.q1.Copy(m.q1s[i])
			m.q1s[i].ScaleD(m.rates[i])
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
		for i := 0; i < m.ncat; i++ {
			m.q2.Copy(m.q2s[i])
			m.q2s[i].ScaleD(m.rates[i])
		}
		m.q2done = true
	}

	m.updateProportions()
	m.expAllBr = false
}

func (m *BranchSiteGamma) updateRates() {
	m.rates = paml.DiscreteGamma(m.alpha, m.alpha, m.ncat, false, nil, m.rates)
	m.q0done = false
	m.q1done = false
	m.q2done = false
	m.ratesdone = true
}

func (m *BranchSiteGamma) Likelihood() float64 {
	if !m.ratesdone {
		m.updateRates()
	}
	if !m.q0done || !m.q1done || !m.q2done {
		m.updateMatrices()
	}
	if !m.propdone {
		m.updateProportions()
	}
	return m.BaseModel.Likelihood()
}

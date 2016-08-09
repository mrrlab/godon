package cmodel

import (
	"math"
	"math/rand"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/tree"
)

// BranchSiteC is a special case of the branch-site model. p2 is
// computed based on omega2 value: p2 = omega2 / 1000.
type BranchSiteC struct {
	*BaseModel
	q0, q1, q2             *codon.EMatrix
	kappa                  float64
	omega0, omega2         float64
	p0prop                 float64
	q0done, q1done, q2done bool
	propdone               bool
}

// NewBranchSiteC creates a new BranchSiteC model.
func NewBranchSiteC(cali codon.CodonSequences, t *tree.Tree, cf codon.CodonFrequency) (m *BranchSiteC) {
	m = &BranchSiteC{
		q0: &codon.EMatrix{CF: cf},
		q1: &codon.EMatrix{CF: cf},
		q2: &codon.EMatrix{CF: cf},
	}

	m.BaseModel = NewBaseModel(cali, t, cf, m)

	m.setupParameters()
	m.SetBranchMatrices()
	m.SetDefaults()

	return
}

// GetNClass returns number of site classes.
func (m *BranchSiteC) GetNClass() int {
	return 4
}

// Copy makes a copy of the model preserving the model parameter
// values.
func (m *BranchSiteC) Copy() optimize.Optimizable {
	newM := &BranchSiteC{
		BaseModel: m.BaseModel.Copy(),
		q0:        &codon.EMatrix{CF: m.cf},
		q1:        &codon.EMatrix{CF: m.cf},
		q2:        &codon.EMatrix{CF: m.cf},
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

// addParameters adds all the model parameters to the parameter
// storage.
func (m *BranchSiteC) addParameters(fpg optimize.FloatParameterGenerator) {
	kappa := fpg(&m.kappa, "kappa")
	kappa.SetOnChange(func() {
		m.q0done = false
		m.q1done = false
		m.q2done = false
	})
	kappa.SetPriorFunc(optimize.UniformPrior(0, 20, false, true))
	kappa.SetProposalFunc(optimize.NormalProposal(0.01))
	kappa.SetMin(1e-2)
	kappa.SetMax(20)
	m.parameters.Append(kappa)

	omega0 := fpg(&m.omega0, "omega0")
	omega0.SetOnChange(func() {
		m.q0done = false
	})
	omega0.SetPriorFunc(optimize.GammaPrior(1, 2, false))
	omega0.SetProposalFunc(optimize.NormalProposal(0.01))
	omega0.SetMin(1e-4)
	omega0.SetMax(1)
	m.parameters.Append(omega0)

	omega2 := fpg(&m.omega2, "omega2")
	omega2.SetOnChange(func() {
		m.q2done = false
	})
	omega2.SetPriorFunc(optimize.GammaPrior(1, 2, false))
	omega2.SetProposalFunc(optimize.NormalProposal(0.01))
	omega2.SetMin(1)
	omega2.SetMax(1000)
	m.parameters.Append(omega2)

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

// SetParameters sets the model parameter values.
func (m *BranchSiteC) SetParameters(kappa float64, omega0, omega2 float64, p0prop float64) {
	m.kappa = kappa
	m.omega0 = omega0

	m.omega2 = omega2

	m.p0prop = p0prop

	m.q0done, m.q1done, m.q2done = false, false, false
}

// GetParameters returns the model parameter values.
func (m *BranchSiteC) GetParameters() (kappa float64, omega0, omega2 float64, p0 float64) {
	return m.kappa, m.omega0, m.omega2, m.p0prop
}

// SetDefaults sets the default initial parameter values.
func (m *BranchSiteC) SetDefaults() {
	kappa := 1e-2 + rand.Float64()*10
	omega0 := 0.2 + 0.1*rand.Float64()
	omega2 := 3.1 + rand.Float64()
	x0 := 1.0 + 0.5*rand.Float64()
	x1 := 0.2 * rand.Float64()
	p0 := math.Exp(x0) / (1 + math.Exp(x0) + math.Exp(x1))
	m.SetParameters(kappa, omega0, omega2, p0)
}

// SetBranchMatrices set matrices for all the branches.
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

// updateProportions updates proportions if model parameters are
// changing.
func (m *BranchSiteC) updateProportions() {
	p2 := m.omega2 / 1000
	m.prop[0][0] = m.p0prop * (1 - p2)
	m.prop[0][1] = (1 - m.p0prop) * (1 - p2)
	m.prop[0][2] = m.p0prop * p2
	m.prop[0][3] = (1 - m.p0prop) * p2

	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		if node.Class == 0 {
			m.scale[node.Id] = (m.prop[0][0]+m.prop[0][2])*m.q0.Scale + (m.prop[0][1]+m.prop[0][3])*m.q1.Scale
		} else {
			m.scale[node.Id] = m.prop[0][0]*m.q0.Scale + m.prop[0][1]*m.q1.Scale + (m.prop[0][2]+m.prop[0][3])*m.q2.Scale
		}
	}
	m.propdone = true
	m.expAllBr = false
}

// updateMatrices updates matrices if model parameters are changing.
func (m *BranchSiteC) updateMatrices() {
	if !m.q0done {
		Q0, s0 := codon.CreateTransitionMatrix(m.cf, m.kappa, m.omega0, m.q0.Q)
		m.q0.Set(Q0, s0)
		err := m.q0.Eigen()
		if err != nil {
			panic("error eigen q0")
		}
		m.q0done = true
	}

	if !m.q1done {
		Q1, s1 := codon.CreateTransitionMatrix(m.cf, m.kappa, 1, m.q1.Q)
		m.q1.Set(Q1, s1)
		err := m.q1.Eigen()
		if err != nil {
			panic("error eigen q1")
		}
		m.q1done = true
	}

	if !m.q2done {
		Q2, s2 := codon.CreateTransitionMatrix(m.cf, m.kappa, m.omega2, m.q2.Q)
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

// Likelihood computes likelihood.
func (m *BranchSiteC) Likelihood() float64 {
	if !m.q0done || !m.q1done || !m.q2done {
		m.updateMatrices()
	}
	if !m.propdone {
		m.updateProportions()
	}
	return m.BaseModel.Likelihood()
}

package cmodel

import (
	"fmt"
	"math"
	"math/rand"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/paml"
)

// BranchSiteGammaERates is an implementation of the branch-site model with
// gamma rates variation with exlicit rate values.
type BranchSiteGammaERates struct {
	*BaseModel
	q0s            []*codon.EMatrix
	q1s            []*codon.EMatrix
	q2s            []*codon.EMatrix
	kappa          float64
	omega0, omega2 float64
	p01sum, p0prop float64
	// site gamma alpha parameter
	alphas float64
	gammas []float64
	// codon gamma alpha parameter
	alphac float64
	gammac []float64
	// this should be discrete, but optimizer supports only float
	// so class[i] = int(cs_rates[i] * NClass / 4)
	csRates []float64
	// temporary array
	tmp                    []float64
	ncatcg                 int
	ncatsg                 int
	fixw2                  bool
	q0done, q1done, q2done bool
	allpropdone            bool
	propdone               []bool
	gammacdone             bool
	gammasdone             bool
	csrdone                bool
}

// NewBranchSiteGammaERates creates a new BranchSiteGammaERates model.
func NewBranchSiteGammaERates(data *Data, fixw2 bool, ncatsg, ncatcg int) (m *BranchSiteGammaERates) {
	scat := ncatsg * ncatsg * ncatsg

	m = &BranchSiteGammaERates{
		fixw2:    fixw2,
		ncatsg:   ncatsg,
		ncatcg:   ncatcg,
		q0s:      make([]*codon.EMatrix, scat*ncatcg),
		q1s:      make([]*codon.EMatrix, scat*ncatcg),
		q2s:      make([]*codon.EMatrix, scat*ncatcg),
		gammas:   make([]float64, ncatsg),
		gammac:   make([]float64, ncatcg),
		csRates:  make([]float64, data.cSeqs.Length()),
		propdone: make([]bool, data.cSeqs.Length()),
		tmp:      make([]float64, maxInt(ncatcg, ncatsg, 3)),
	}
	m.BaseModel = NewBaseModel(data, m)

	// We need to create category for every site and every codon
	for i := 0; i < data.cSeqs.Length(); i++ {
		m.prop[i] = make([]float64, m.GetNClass())
	}

	for i := 0; i < scat*ncatcg; i++ {
		m.q0s[i] = &codon.EMatrix{CF: data.cFreq}
		m.q1s[i] = &codon.EMatrix{CF: data.cFreq}
		m.q2s[i] = &codon.EMatrix{CF: data.cFreq}
	}

	m.setupParameters()
	m.setBranchMatrices()
	m.SetDefaults()

	return

}

// GetNClass returns number of site classes.
func (m *BranchSiteGammaERates) GetNClass() int {
	scat := m.ncatsg * m.ncatsg * m.ncatsg

	return 4 * scat * m.ncatcg
}

// Copy makes a copy of the model preserving the model parameter
// values.
func (m *BranchSiteGammaERates) Copy() optimize.Optimizable {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	newM := &BranchSiteGammaERates{
		BaseModel: m.BaseModel.Copy(),
		q0s:       make([]*codon.EMatrix, scat*m.ncatcg),
		q1s:       make([]*codon.EMatrix, scat*m.ncatcg),
		q2s:       make([]*codon.EMatrix, scat*m.ncatcg),
		kappa:     m.kappa,
		omega0:    m.omega0,
		omega2:    m.omega2,
		p01sum:    m.p01sum,
		p0prop:    m.p0prop,
		fixw2:     m.fixw2,

		ncatsg: m.ncatsg,
		alphas: m.alphas,
		gammas: make([]float64, m.ncatsg),

		ncatcg: m.ncatcg,
		alphac: m.alphac,
		gammac: make([]float64, m.ncatcg),

		tmp: make([]float64, maxInt(m.ncatcg, m.ncatsg, 3)),
	}
	newM.BaseModel.Model = newM

	copy(newM.csRates, m.csRates)

	// basemodel assumes that all props are identical,
	// we need to copy them explicitly
	for i := 0; i < m.data.cSeqs.Length(); i++ {
		newM.prop[i] = make([]float64, m.GetNClass())
		copy(newM.prop[i], m.prop[i])
	}

	for i := 0; i < scat*m.ncatcg; i++ {
		newM.q0s[i] = &codon.EMatrix{CF: m.data.cFreq}
		newM.q1s[i] = &codon.EMatrix{CF: m.data.cFreq}
		newM.q2s[i] = &codon.EMatrix{CF: m.data.cFreq}
	}

	newM.setupParameters()
	newM.setBranchMatrices()
	return newM
}

// addParameters adds all the model parameters to the parameter
// storage.
func (m *BranchSiteGammaERates) addParameters(fpg optimize.FloatParameterGenerator) {
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

	if !m.fixw2 {
		omega2 := fpg(&m.omega2, "omega2")
		omega2.SetOnChange(func() {
			m.q2done = false
		})
		omega2.SetPriorFunc(optimize.GammaPrior(1, 2, false))
		omega2.SetProposalFunc(optimize.NormalProposal(0.01))
		omega2.SetMin(1)
		omega2.SetMax(1000)
		m.parameters.Append(omega2)
	}

	p01sum := fpg(&m.p01sum, "p01sum")
	p01sum.SetOnChange(func() {
		m.allpropdone = false
	})
	p01sum.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
	p01sum.SetMin(1e-12)
	p01sum.SetMax(1)
	p01sum.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(p01sum)

	p0prop := fpg(&m.p0prop, "p0prop")
	p0prop.SetOnChange(func() {
		m.allpropdone = false
	})
	p0prop.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
	p0prop.SetMin(0)
	p0prop.SetMax(1)
	p0prop.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(p0prop)

	if m.ncatsg > 1 {
		alphas := fpg(&m.alphas, "alphas")
		alphas.SetOnChange(func() {
			m.gammasdone = false
		})
		alphas.SetPriorFunc(optimize.GammaPrior(1, 2, false))
		alphas.SetMin(1e-2)
		alphas.SetMax(1000)
		alphas.SetProposalFunc(optimize.NormalProposal(0.01))
		m.parameters.Append(alphas)
	}

	if m.ncatcg > 1 {
		alphac := fpg(&m.alphac, "alphac")
		alphac.SetOnChange(func() {
			m.gammacdone = false
		})
		alphac.SetPriorFunc(optimize.GammaPrior(1, 2, false))
		alphac.SetMin(1e-2)
		alphac.SetMax(1000)
		alphac.SetProposalFunc(optimize.NormalProposal(0.01))
		m.parameters.Append(alphac)

	}

	if m.ncatcg > 1 || m.ncatsg > 1 {
		// every codon has a rate parameter
		for i := 0; i < m.data.cSeqs.Length(); i++ {
			nm := fmt.Sprintf("rate_%04d", i+1)
			rate := optimize.NewDiscreteParameter(&m.csRates[i], nm, m.GetNClass()/4)
			// make a copy of i for the closure
			pos := i
			rate.SetOnChange(func() {
				m.propdone[pos] = false
				m.csrdone = false
			})

			m.parameters.Append(rate)
		}
	}

}

// SetParameters sets the model parameter values.
func (m *BranchSiteGammaERates) SetParameters(kappa float64, omega0, omega2 float64, p0, p1 float64, alphas, alphac float64, csRates []float64) {
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
	m.alphas = alphas
	m.alphac = alphac
	if len(csRates) != len(m.csRates) {
		panic("Incorrect number of rates")
	}
	copy(m.csRates, csRates)
}

// GetParameters returns the model parameter values.
func (m *BranchSiteGammaERates) GetParameters() (kappa float64, omega0, omega2 float64, p0, p1 float64, alphas, alphac float64, csRates []float64) {
	csRates = make([]float64, len(m.csRates))
	copy(csRates, m.csRates)
	return m.kappa, m.omega0, m.omega2, m.p01sum * m.p0prop, m.p01sum * (1 - m.p0prop), m.alphas, m.alphac, csRates
}

// SetDefaults sets the default initial parameter values.
func (m *BranchSiteGammaERates) SetDefaults() {
	// these parameters mostly come from codeml
	kappa := 1e-2 + rand.Float64()*10
	omega0 := 0.2 + 0.1*rand.Float64()
	omega2 := 3.1 + rand.Float64()
	x0 := 1.0 + 0.5*rand.Float64()
	x1 := 0.2 * rand.Float64()
	p0 := math.Exp(x0) / (1 + math.Exp(x0) + math.Exp(x1))
	p1 := math.Exp(x1) / (1 + math.Exp(x0) + math.Exp(x1))
	alphas := 0.5 + rand.Float64()*3
	alphac := 0.5 + rand.Float64()*3

	csRates := make([]float64, m.data.cSeqs.Length())
	for i := range csRates {
		csRates[i] = float64(rand.Intn(m.GetNClass() / 4))
	}

	m.SetParameters(kappa, omega0, omega2, p0, p1, alphas, alphac, csRates)
}

// Organization of the class categories.
// cat 1, internal gamma cat 1, external gamma cat 1
// cat 1, internal gamma cat 1, external gamma cat 2
// ...
// cat 1, internal gamma cat 1, external gamma ncatcg
// ...
// cat 1, internal gamma cat ncatsg^3, external gamma cat 1
// ...
// cat 1, internal gamma cat ncatsg^3, external gamma ncatcg
// cat 2, internal gamma cat 1, external gamma cat 1
// ...
// cat 2, internal gamma cat ncatsg^3, external gamma cat ncatcg
// ...
// ...
// cat 4, internal gamma cat 1, external gamma cat 1
// ...
// cat 4, internal gamma cat ncatsg^3, external gamma cat ncatcg
// (total: 4 * ncatsg^3 * ncatcg)

// setBranchMatrices set matrices for all the branches.
func (m *BranchSiteGammaERates) setBranchMatrices() {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	for _, node := range m.data.Tree.NodeIDArray() {
		if node == nil {
			continue
		}
		for i := 0; i < m.ncatcg; i++ {
			for j := 0; j < scat; j++ {
				m.qs[i+j*m.ncatcg+0*scat*m.ncatcg][node.ID] = m.q0s[i+j*m.ncatcg]
				m.qs[i+j*m.ncatcg+1*scat*m.ncatcg][node.ID] = m.q1s[i+j*m.ncatcg]
				if node.Class == 1 {
					m.qs[i+j*m.ncatcg+2*scat*m.ncatcg][node.ID] = m.q2s[i+j*m.ncatcg]
					m.qs[i+j*m.ncatcg+3*scat*m.ncatcg][node.ID] = m.q2s[i+j*m.ncatcg]
				} else {
					m.qs[i+j*m.ncatcg+2*scat*m.ncatcg][node.ID] = m.q0s[i+j*m.ncatcg]
					m.qs[i+j*m.ncatcg+3*scat*m.ncatcg][node.ID] = m.q1s[i+j*m.ncatcg]
				}
			}
		}
	}
}

// updateProportions updates proportions if model parameters are
// changing.
func (m *BranchSiteGammaERates) updateProportions() {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	p0 := m.p0prop * m.p01sum
	p1 := m.p01sum - p0
	bothcat := m.ncatcg * scat

	for pos := 0; pos < m.data.cSeqs.Length(); pos++ {
		if !m.allpropdone || !m.propdone[pos] {
			class := int(m.csRates[pos])
			for i := 0; i < m.GetNClass()/4; i++ {
				if class == i {
					m.prop[pos][i+bothcat*0] = p0
					m.prop[pos][i+bothcat*1] = p1
					m.prop[pos][i+bothcat*2] = (1 - p0 - p1) * p0 / (p0 + p1)
					m.prop[pos][i+bothcat*3] = (1 - p0 - p1) * p1 / (p0 + p1)

				} else {
					m.prop[pos][i+bothcat*0] = 0
					m.prop[pos][i+bothcat*1] = 0
					m.prop[pos][i+bothcat*2] = 0
					m.prop[pos][i+bothcat*3] = 0
				}
			}
			m.prunPos[pos] = false
		}
	}

	for _, node := range m.data.Tree.NodeIDArray() {
		if node == nil {
			continue
		}
		scale := 0.0
		for i := 0; i < m.data.cSeqs.Length(); i++ {
			for j := 0; j < m.GetNClass(); j++ {
				scale += m.prop[i][j] * m.qs[j][node.ID].Scale
			}
		}
		m.scale[node.ID] = scale / float64(m.data.cSeqs.Length())
	}
	if !m.allpropdone {
		m.prunAllPos = false
	}
	m.allpropdone = true
	m.expAllBr = false
}

// fillMatrices sets Q-matrices for all the rates.
func (m *BranchSiteGammaERates) fillMatrices(omega float64, dest []*codon.EMatrix) {
	for c1 := 0; c1 < m.ncatsg; c1++ {
		m.tmp[0] = m.gammas[c1]
		for c2 := 0; c2 < m.ncatsg; c2++ {
			m.tmp[1] = m.gammas[c2]
			for c3 := 0; c3 < m.ncatsg; c3++ {
				m.tmp[2] = m.gammas[c3]

				e := &codon.EMatrix{CF: m.data.cFreq}
				Q, s := codon.CreateRateTransitionMatrix(m.data.cFreq, m.kappa, omega, m.tmp, e.Q)
				e.Set(Q, s)
				err := e.Eigen()
				if err != nil {
					panic("error finding eigen")
				}

				for ecl, rate := range m.gammac {
					catid := (((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl

					e.Copy(dest[catid])
					dest[catid].ScaleD(rate)
				}
			}
		}

	}
}

// updateMatrices updates matrices if model parameters are changing.
func (m *BranchSiteGammaERates) updateMatrices() {
	if !m.q0done {
		m.fillMatrices(m.omega0, m.q0s)
		m.q0done = true
	}

	if !m.q1done {
		m.fillMatrices(1, m.q1s)
		m.q1done = true
	}

	if !m.q2done {
		m.fillMatrices(m.omega2, m.q2s)
		m.q2done = true
	}

	m.updateProportions()
	m.expAllBr = false
}

// Likelihood computes likelihood.
func (m *BranchSiteGammaERates) Likelihood() float64 {
	if !m.gammasdone {
		if m.ncatsg > 1 {
			m.gammas = paml.DiscreteGamma(m.alphas, m.alphas, m.ncatsg, false, m.tmp, m.gammas)
			m.q0done = false
			m.q1done = false
			m.q2done = false
		} else {
			m.gammas[0] = 1
		}
		m.gammasdone = true
	}
	if !m.gammacdone {
		if m.ncatcg > 1 {
			m.gammac = paml.DiscreteGamma(m.alphac, m.alphac, m.ncatcg, false, m.tmp, m.gammac)
			m.q0done = false
			m.q1done = false
			m.q2done = false
		} else {
			m.gammac[0] = 1
		}
		m.gammacdone = true
	}
	if !m.q0done || !m.q1done || !m.q2done {
		m.updateMatrices()
	}
	if !m.allpropdone || !m.csrdone {
		m.updateProportions()
	}
	return m.BaseModel.Likelihood()
}

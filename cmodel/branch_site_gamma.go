package cmodel

import (
	"math"
	"math/rand"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/paml"
	"bitbucket.org/Davydov/godon/tree"
)

// BranchSite is an implementation of the branch-site model with gamma
// rates variation.
type BranchSiteGamma struct {
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
	// temporary array
	tmp                    []float64
	ncatcg                 int
	ncatsg                 int
	fixw2                  bool
	q0done, q1done, q2done bool
	propdone               bool
	gammacdone             bool
	gammasdone             bool
}

// NewBranchSiteGamma creates a new BranchSiteGamma model.
func NewBranchSiteGamma(cali codon.CodonSequences, t *tree.Tree, cf codon.CodonFrequency, fixw2 bool, ncatsg, ncatcg int) (m *BranchSiteGamma) {
	scat := ncatsg * ncatsg * ncatsg

	m = &BranchSiteGamma{
		fixw2:  fixw2,
		ncatsg: ncatsg,
		ncatcg: ncatcg,
		q0s:    make([]*codon.EMatrix, scat*ncatcg),
		q1s:    make([]*codon.EMatrix, scat*ncatcg),
		q2s:    make([]*codon.EMatrix, scat*ncatcg),
		gammas: make([]float64, ncatsg),
		gammac: make([]float64, ncatcg),
		tmp:    make([]float64, maxInt(ncatcg, ncatsg, 3)),
	}
	m.BaseModel = NewBaseModel(cali, t, cf, m)

	for i := 0; i < scat*ncatcg; i++ {
		m.q0s[i] = &codon.EMatrix{CF: cf}
		m.q1s[i] = &codon.EMatrix{CF: cf}
		m.q2s[i] = &codon.EMatrix{CF: cf}
	}

	m.setupParameters()
	m.setBranchMatrices()
	m.SetDefaults()

	return

}

// GetNClass returns number of site classes.
func (m *BranchSiteGamma) GetNClass() int {
	scat := m.ncatsg * m.ncatsg * m.ncatsg

	return 4 * scat * m.ncatcg
}

// Copy makes a copy of the model preserving the model parameter
// values.
func (m *BranchSiteGamma) Copy() optimize.Optimizable {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	newM := &BranchSiteGamma{
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

	for i := 0; i < scat*m.ncatcg; i++ {
		newM.q0s[i] = &codon.EMatrix{CF: m.cf}
		newM.q1s[i] = &codon.EMatrix{CF: m.cf}
		newM.q2s[i] = &codon.EMatrix{CF: m.cf}
	}

	newM.setupParameters()
	newM.setBranchMatrices()
	return newM
}

// addParameters adds all the model parameters to the parameter
// storage.
func (m *BranchSiteGamma) addParameters(fpg optimize.FloatParameterGenerator) {
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
}

// SetParameters sets the model parameter values.
func (m *BranchSiteGamma) SetParameters(kappa float64, omega0, omega2 float64, p0, p1 float64, alphas, alphac float64) {
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
}

// GetParameters returns the model parameter values.
func (m *BranchSiteGamma) GetParameters() (kappa float64, omega0, omega2 float64, p0, p1 float64, alphas, alphac float64) {
	return m.kappa, m.omega0, m.omega2, m.p01sum * m.p0prop, m.p01sum * (1 - m.p0prop), m.alphas, m.alphac
}

// SetDefaults sets the default initial parameter values.
func (m *BranchSiteGamma) SetDefaults() {
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
	m.SetParameters(kappa, omega0, omega2, p0, p1, alphas, alphac)
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
func (m *BranchSiteGamma) setBranchMatrices() {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		for i := 0; i < m.ncatcg; i++ {
			for j := 0; j < scat; j++ {
				m.qs[i+j*m.ncatcg+0*scat*m.ncatcg][node.Id] = m.q0s[i+j*m.ncatcg]
				m.qs[i+j*m.ncatcg+1*scat*m.ncatcg][node.Id] = m.q1s[i+j*m.ncatcg]
				if node.Class == 1 {
					m.qs[i+j*m.ncatcg+2*scat*m.ncatcg][node.Id] = m.q2s[i+j*m.ncatcg]
					m.qs[i+j*m.ncatcg+3*scat*m.ncatcg][node.Id] = m.q2s[i+j*m.ncatcg]
				} else {
					m.qs[i+j*m.ncatcg+2*scat*m.ncatcg][node.Id] = m.q0s[i+j*m.ncatcg]
					m.qs[i+j*m.ncatcg+3*scat*m.ncatcg][node.Id] = m.q1s[i+j*m.ncatcg]
				}
			}
		}
	}
}

// updateProportions updates proportions if model parameters are
// changing.
func (m *BranchSiteGamma) updateProportions() {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	p0 := m.p0prop * m.p01sum
	p1 := m.p01sum - p0
	bothcat := m.ncatcg * scat
	for i := 0; i < bothcat; i++ {
		m.prop[0][i+bothcat*0] = p0 / float64(bothcat)
		m.prop[0][i+bothcat*1] = p1 / float64(bothcat)
		m.prop[0][i+bothcat*2] = (1 - p0 - p1) * p0 / (p0 + p1) / float64(bothcat)
		m.prop[0][i+bothcat*3] = (1 - p0 - p1) * p1 / (p0 + p1) / float64(bothcat)
	}

	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		scale := 0.0
		for i := 0; i < bothcat*4; i++ {
			scale += m.prop[0][i] * m.qs[i][node.Id].Scale
		}
		m.scale[node.Id] = scale
	}
	m.propdone = true
	m.expAllBr = false
}

// fillMatrices sets Q-matrices for all the rates.
func (m *BranchSiteGamma) fillMatrices(omega float64, dest []*codon.EMatrix) {
	for c1 := 0; c1 < m.ncatsg; c1++ {
		m.tmp[0] = m.gammas[c1]
		for c2 := 0; c2 < m.ncatsg; c2++ {
			m.tmp[1] = m.gammas[c2]
			for c3 := 0; c3 < m.ncatsg; c3++ {
				m.tmp[2] = m.gammas[c3]

				e := &codon.EMatrix{CF: m.cf}
				Q, s := codon.CreateRateTransitionMatrix(m.cf, m.kappa, omega, m.tmp, e.Q)
				e.Set(Q, s)
				err := e.Eigen()
				if err != nil {
					log.Fatal(err)
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
func (m *BranchSiteGamma) updateMatrices() {
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

	m.propdone = false
	m.expAllBr = false
}

// Final prints NEB results (only if with positive selection).
func (m *BranchSiteGamma) Final() {
	// if w2=1, do not perform NEB analysis.
	if m.fixw2 {
		log.Info("No NEB since no positive selection in the model.")
		return
	}
	classes := make(map[int]bool, m.GetNClass())
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	for i := 0; i < m.ncatcg; i++ {
		for j := 0; j < scat; j++ {
			classes[m.ncatcg*scat*2+i+j*m.ncatcg] = true
			classes[m.ncatcg*scat*3+i+j*m.ncatcg] = true
		}
	}

	posterior := m.NEBPosterior(classes)

	log.Notice("NEB analysis")
	m.PrintPosterior(posterior)

	//posterior = m.BEBPosterior()

	log.Notice("BEB analysis")
	m.PrintPosterior(posterior)
}

// Likelihood computes likelihood.
func (m *BranchSiteGamma) Likelihood() float64 {
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
	if !m.propdone {
		m.updateProportions()
	}
	return m.BaseModel.Likelihood()
}

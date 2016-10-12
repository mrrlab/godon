package cmodel

import (
	"math/rand"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/paml"
	"bitbucket.org/Davydov/godon/tree"
)

// M1 is an implementation of M1a & M2a models.
type M2 struct {
	*BaseModel
	q0             []*codon.EMatrix
	q1             []*codon.EMatrix
	q2             []*codon.EMatrix
	p0, p1prop     float64
	omega0, omega2 float64
	kappa          float64
	// site gamma alpha parameter
	alphas float64
	gammas []float64
	// codon gamma alpha parameter
	alphac float64
	gammac []float64
	// only allow omega2 if addw is true
	addw       bool
	ncatsg     int
	ncatcg     int
	q0done     bool
	q1done     bool
	q2done     bool
	propdone   bool
	gammasdone bool
	gammacdone bool
	// temporary array
	tmp []float64
}

// NewM8 creates a new M8 model.
func NewM2(cali codon.CodonSequences, t *tree.Tree, cf codon.CodonFrequency, addw bool, ncatsg, ncatcg int) (m *M2) {
	// n site gamma categories, ncatb * n^3 matrices
	gcat := ncatsg * ncatsg * ncatsg
	m = &M2{
		addw:   addw,
		ncatsg: ncatsg,
		ncatcg: ncatcg,
		q0:     make([]*codon.EMatrix, gcat*ncatcg),
		q1:     make([]*codon.EMatrix, gcat*ncatcg),
		q2:     make([]*codon.EMatrix, gcat*ncatcg),
		gammas: make([]float64, ncatsg),
		gammac: make([]float64, ncatcg),
		tmp:    make([]float64, maxInt(ncatsg, ncatcg, 3)),
	}

	for i := 0; i < gcat*ncatcg; i++ {
		m.q0[i] = &codon.EMatrix{CF: cf}
		m.q1[i] = &codon.EMatrix{CF: cf}
		if addw {
			m.q2[i] = &codon.EMatrix{CF: cf}
		}
	}

	m.BaseModel = NewBaseModel(cali, t, cf, m)

	m.setupParameters()
	m.setBranchMatrices()
	m.SetDefaults()
	return
}

// GetNClass returns number of site classes.
func (m *M2) GetNClass() int {
	gcat := m.ncatsg * m.ncatsg * m.ncatsg
	if m.addw {
		return gcat * m.ncatcg * 3
	}
	return gcat * m.ncatcg * 2
}

// Copy makes a copy of the model preserving the model parameter
// values.
func (m *M2) Copy() optimize.Optimizable {
	// n inner gamma categories, ncatb * n^3 matrices
	gcat := m.ncatsg * m.ncatsg * m.ncatsg
	newM := &M2{
		BaseModel: m.BaseModel.Copy(),
		q0:        make([]*codon.EMatrix, gcat*m.ncatcg),
		q1:        make([]*codon.EMatrix, gcat*m.ncatcg),
		q2:        make([]*codon.EMatrix, gcat*m.ncatcg),
		tmp:       make([]float64, maxInt(m.ncatsg, m.ncatcg, 3)),
		ncatsg:    m.ncatsg,
		ncatcg:    m.ncatcg,
		gammas:    make([]float64, m.ncatsg),
		gammac:    make([]float64, m.ncatcg),
		addw:      m.addw,
		p0:        m.p0,
		p1prop:    m.p1prop,
		omega0:    m.omega0,
		omega2:    m.omega2,
		kappa:     m.kappa,
		alphas:    m.alphas,
		alphac:    m.alphac,
	}

	for i := 0; i < gcat*m.ncatcg; i++ {
		newM.q0[i] = &codon.EMatrix{CF: m.cf}
		newM.q1[i] = &codon.EMatrix{CF: m.cf}
		if m.addw {
			newM.q2[i] = &codon.EMatrix{CF: m.cf}
		}
	}

	newM.BaseModel.Model = newM
	newM.setupParameters()
	newM.setBranchMatrices()
	return newM
}

// addParameters adds all the model parameters to the parameter
// storage.
func (m *M2) addParameters(fpg optimize.FloatParameterGenerator) {
	p0 := fpg(&m.p0, "p0")
	p0.SetOnChange(func() {
		m.propdone = false
	})
	p0.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
	p0.SetMin(0)
	p0.SetMax(1)
	p0.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(p0)

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
	omega0.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
	omega0.SetProposalFunc(optimize.NormalProposal(0.01))
	omega0.SetMin(0)
	omega0.SetMax(1)
	m.parameters.Append(omega0)

	if m.addw {
		p1prop := fpg(&m.p1prop, "p1prop")
		p1prop.SetOnChange(func() {
			m.propdone = false
		})
		p1prop.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
		p1prop.SetMin(0)
		p1prop.SetMax(1)
		p1prop.SetProposalFunc(optimize.NormalProposal(0.01))
		m.parameters.Append(p1prop)

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

// GetParameters returns the model parameter values.
func (m *M2) GetParameters() (p0, p1prop, omega0, omega2, kappa, alphas, alphac float64) {
	return m.p0, m.p1prop, m.omega0, m.omega2, m.kappa, m.alphas, m.alphac
}

// SetParameters sets the model parameter values.
func (m *M2) SetParameters(p0, p1prop, omega0, omega2, kappa, alphas, alphac float64) {
	if m.addw {
		m.p1prop = p1prop
	} else {
		m.p1prop = 1
	}
	m.p0 = p0
	m.kappa = kappa
	m.omega0 = omega0
	m.omega2 = omega2
	m.alphas = alphas
	m.alphac = alphac
	m.gammasdone = false
	m.gammacdone = false
	m.q0done = false
	m.q1done = false
	m.q2done = false
}

// SetDefaults sets the default initial parameter values.
func (m *M2) SetDefaults() {
	p0 := 0.69 + rand.Float64()*0.3
	p1prop := 0.69 + rand.Float64()*0.3
	kappa := 1e-2 + rand.Float64()*10
	omega0 := 0.2 + 0.1*rand.Float64()
	omega2 := 1.001 + rand.Float64()*10
	alphas := 1e-3 + rand.Float64()*10
	alphac := 1e-3 + rand.Float64()*10
	m.SetParameters(p0, p1prop, omega0, omega2, kappa, alphas, alphac)
}

// Organization of the class categories.
// omega0, internal gamma cat 1, external gamma cat 1
// omega0, internal gamma cat 1, external gamma cat 2
// ...
// omega0, internal gamma cat 1, external gamma ncatcg
// ...
// omega0, internal gamma cat ncatsg^3, external gamma cat 1
// ...
// omega0, internal gamma cat ncatsg^3, external gamma ncatcg
// omega1, internal gamma cat 1, external gamma cat 1
// ...
// omega1, internal gamma cat ncatsg^3, external gamma cat ncatcg
// (total: 2 * ncatsg^3 * ncatcg)
// if m.addw == true [
//   omega2, internal gamma cat 1, external gamma cat 1
//   ...
//   omega2, internal gamma cat ncatsg^3, external gamma cat ncatcg
// ]
// (total: 3 * ncatsg^3 * ncatcg

// setBranchMatrices set matrices for all the branches.
func (m *M2) setBranchMatrices() {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		for i := 0; i < m.ncatcg; i++ {
			for j := 0; j < scat; j++ {
				m.qs[i+j*m.ncatcg+0*scat*m.ncatcg][node.Id] = m.q0[i+j*m.ncatcg]
				m.qs[i+j*m.ncatcg+1*scat*m.ncatcg][node.Id] = m.q1[i+j*m.ncatcg]
				if m.addw {
					m.qs[i+j*m.ncatcg+2*scat*m.ncatcg][node.Id] = m.q2[i+j*m.ncatcg]
				}
			}
		}
	}
}

// fillMatrices sets Q-matrices for all the rates.
func (m *M2) fillMatrices(omega float64, dest []*codon.EMatrix) {
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
func (m *M2) updateMatrices() {
	if !m.q0done {
		m.fillMatrices(m.omega0, m.q0)
		m.q0done = true
	}

	if !m.q1done {
		m.fillMatrices(1, m.q1)
		m.q1done = true
	}

	if m.addw && !m.q2done {
		m.fillMatrices(m.omega2, m.q2)
		m.q2done = true
	}

	m.propdone = false
}

// Final prints NEB results (only if with positive selection).
func (m *M2) Final() {
	// if w2=1, do not perform NEB analysis.
	if !m.addw {
		log.Info("No NEB since no positive selection in the model.")
		return
	}
	classes := make(map[int]bool, m.GetNClass())

	gcat := m.ncatsg * m.ncatsg * m.ncatsg

	for c1 := 0; c1 < m.ncatsg; c1++ {
		for c2 := 0; c2 < m.ncatsg; c2++ {
			for c3 := 0; c3 < m.ncatsg; c3++ {
				for ecl := range m.gammac {
					catid := (((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl

					class := 2*gcat*m.ncatcg + catid
					classes[class] = true
				}
			}
		}
	}

	posterior := m.NEBPosterior(classes)

	m.PrintPosterior(posterior)
}

// updateProportions updates proportions if model parameters are
// changing.
func (m *M2) updateProportions() {
	scat := m.ncatsg * m.ncatsg * m.ncatsg

	p0 := m.p0
	p1 := (1 - m.p0) * m.p1prop
	p2 := 1 - p0 - p1

	bothcat := m.ncatcg * scat
	for i := 0; i < bothcat; i++ {
		m.prop[0][i+bothcat*0] = p0 / float64(bothcat)
		m.prop[0][i+bothcat*1] = p1 / float64(bothcat)
		if m.addw {
			m.prop[0][i+bothcat*2] = p2 / float64(bothcat)
		}
	}

	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		scale := 0.0
		for i := range m.prop[0] {
			scale += m.prop[0][i] * m.qs[i][node.Id].Scale
		}
		m.scale[node.Id] = scale
	}
	m.propdone = true
	m.expAllBr = false
}

// Likelihood computes likelihood.
func (m *M2) Likelihood() float64 {
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

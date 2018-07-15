package cmodel

import (
	"math"
	"math/rand"
	"time"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/dist"
	"bitbucket.org/Davydov/godon/optimize"
)

// M8 is an implementation of M8 model.
type M8 struct {
	*BaseModel
	qb           []*codon.EMatrix
	q0           []*codon.EMatrix
	p0           float64
	p, q         float64
	omega, kappa float64
	// site gamma alpha parameter
	alphas     float64
	gammas     []float64
	gammasprop []float64
	// parametrization inspired by Scheffler 2006
	// proportions
	ps1s, ps2s float64
	// rates, rs2s is inverted
	rs1s, rs2s float64
	// codon gamma alpha parameter
	alphac     float64
	gammac     []float64
	gammacprop []float64
	// parametrization inspired by Scheffler 2006
	// proportions
	ps1c, ps2c float64
	// rates, rs2c is inverted
	rs1c, rs2c float64
	// omega values for beta distribution
	omegab []float64
	// temporary array for beta computations
	tmp []float64
	// only allow extra omega if addw is true
	addw bool
	// fix omega=1
	fixw         bool
	proportional bool
	ncatb        int
	ncatsg       int
	ncatcg       int
	q0done       bool
	qbdone       bool
	propdone     bool
	gammasdone   bool
	gammacdone   bool
	summary      m8Summary
}

// m8Summary stores summary information.
type m8Summary struct {
	SitePosteriorNEB []float64 `json:"sitePosteriorNEB,omitempty"`
	CodonGammaRates  []float64 `json:"codonGammaRates,omitempty"`
	SiteGammaRates   []float64 `json:"siteGammaRates,omitempty"`
	CodonOmega       []float64 `json:"codonOmega,omitempty"`
	PosteriorTime    float64   `json:"posteriorTime,omitempty"`
}

// Empty returns true if there's no data in the structure.
func (s m8Summary) Empty() bool {
	if s.SitePosteriorNEB == nil && s.CodonGammaRates == nil && s.CodonOmega == nil {
		return true
	}
	return false
}

// NewM8 creates a new M8 model.
func NewM8(data *Data, addw, fixw bool, ncatb, ncatsg, ncatcg int, proportional bool) (m *M8) {
	// n site gamma categories, ncatb * n^3 matrices
	gcat := ncatsg * ncatsg * ncatsg
	if ncatb < 2 {
		panic("M8 requires at least two categories")
	}
	m = &M8{
		addw:         addw,
		fixw:         fixw,
		proportional: proportional,
		ncatb:        ncatb,
		ncatsg:       ncatsg,
		ncatcg:       ncatcg,
		qb:           make([]*codon.EMatrix, ncatb*gcat*ncatcg),
		q0:           make([]*codon.EMatrix, gcat*ncatcg),
		gammas:       make([]float64, ncatsg),
		gammasprop:   make([]float64, ncatsg),
		gammac:       make([]float64, ncatcg),
		gammacprop:   make([]float64, ncatcg),
		tmp:          make([]float64, maxInt(ncatb, ncatsg, ncatcg, 3)),
	}

	for i := 0; i < gcat*ncatcg; i++ {
		m.q0[i] = codon.NewEMatrix(data.cFreq)
	}
	for i := 0; i < gcat*ncatcg*ncatb; i++ {
		m.qb[i] = codon.NewEMatrix(data.cFreq)
	}

	m.BaseModel = NewBaseModel(data, m)

	if !m.proportional {
		for i := range m.gammacprop {
			m.gammacprop[i] = 1 / float64(len(m.gammacprop))
		}
		for i := range m.gammasprop {
			m.gammasprop[i] = 1 / float64(len(m.gammasprop))
		}
	}

	m.setupParameters()
	m.SetDefaults()
	return
}

// GetNClass returns number of site classes.
func (m *M8) GetNClass() int {
	gcat := m.ncatsg * m.ncatsg * m.ncatsg
	if m.addw {
		return (gcat * (m.ncatb + 1)) * m.ncatcg
	}
	return gcat * m.ncatb * m.ncatcg
}

// Copy makes a copy of the model preserving the model parameter
// values.
func (m *M8) Copy() optimize.Optimizable {
	// n inner gamma categories, ncatb * n^3 matrices
	gcat := m.ncatsg * m.ncatsg * m.ncatsg
	newM := &M8{
		BaseModel:    m.BaseModel.Copy(),
		qb:           make([]*codon.EMatrix, m.ncatb*gcat*m.ncatcg),
		q0:           make([]*codon.EMatrix, gcat*m.ncatcg),
		tmp:          make([]float64, maxInt(m.ncatb, m.ncatsg, m.ncatcg, 3)),
		ncatb:        m.ncatb,
		proportional: m.proportional,
		ncatsg:       m.ncatsg,
		ncatcg:       m.ncatcg,
		gammas:       make([]float64, m.ncatsg),
		gammasprop:   make([]float64, m.ncatsg),
		gammac:       make([]float64, m.ncatcg),
		gammacprop:   make([]float64, m.ncatcg),
		addw:         m.addw,
		fixw:         m.fixw,
		p0:           m.p0,
		p:            m.p,
		q:            m.q,
		omega:        m.omega,
		kappa:        m.kappa,
		alphas:       m.alphas,
		alphac:       m.alphac,
	}

	for i := 0; i < gcat*m.ncatcg; i++ {
		newM.q0[i] = codon.NewEMatrix(m.data.cFreq)
	}
	for i := 0; i < gcat*m.ncatb*m.ncatcg; i++ {
		newM.qb[i] = codon.NewEMatrix(m.data.cFreq)
	}

	if !m.proportional {
		for i := range m.gammacprop {
			m.gammacprop[i] = 1 / float64(len(m.gammacprop))
		}
		for i := range m.gammasprop {
			m.gammasprop[i] = 1 / float64(len(m.gammasprop))
		}
	}

	newM.BaseModel.model = newM
	newM.setupParameters()
	return newM
}

// addParameters adds all the model parameters to the parameter
// storage.
func (m *M8) addParameters(fpg optimize.FloatParameterGenerator) {
	if m.addw {
		p0 := fpg(&m.p0, "p0")
		p0.SetOnChange(func() {
			m.propdone = false
		})
		p0.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
		p0.SetMin(0)
		p0.SetMax(1)
		p0.SetProposalFunc(optimize.NormalProposal(0.01))
		m.parameters.Append(p0)
	}

	p := fpg(&m.p, "p")
	p.SetOnChange(func() {
		m.qbdone = false
	})
	p.SetPriorFunc(optimize.ExponentialPrior(1, false))
	p.SetMin(0.005)
	p.SetMax(100)
	p.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(p)

	q := fpg(&m.q, "q")
	q.SetOnChange(func() {
		m.qbdone = false
	})
	q.SetPriorFunc(optimize.ExponentialPrior(1, false))
	q.SetMin(0.005)
	q.SetMax(100)
	q.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(q)

	kappa := fpg(&m.kappa, "kappa")
	kappa.SetOnChange(func() {
		m.q0done = false
		m.qbdone = false
	})
	kappa.SetPriorFunc(optimize.UniformPrior(0, 20, false, true))
	kappa.SetProposalFunc(optimize.NormalProposal(0.01))
	kappa.SetMin(1e-2)
	kappa.SetMax(100)
	m.parameters.Append(kappa)

	if m.addw && !m.fixw {
		omega := fpg(&m.omega, "omega")
		omega.SetOnChange(func() {
			m.q0done = false
		})
		omega.SetPriorFunc(optimize.GammaPrior(1, 2, false))
		omega.SetProposalFunc(optimize.NormalProposal(0.01))
		omega.SetMin(1)
		omega.SetMax(1000)
		m.parameters.Append(omega)
	}

	if m.ncatsg > 1 {
		if !m.proportional {
			alphas := fpg(&m.alphas, "alphas")
			alphas.SetOnChange(func() {
				m.gammasdone = false
			})
			alphas.SetPriorFunc(optimize.GammaPrior(1, 2, false))
			alphas.SetMin(1e-2)
			alphas.SetMax(1000)
			alphas.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(alphas)
		} else {
			ps1s := fpg(&m.ps1s, "ps1s")
			ps1s.SetOnChange(func() {
				m.gammasdone = false
			})
			ps1s.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			ps1s.SetMin(1e-5)
			ps1s.SetMax(1 - 1e-5)
			ps1s.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(ps1s)
			ps2s := fpg(&m.ps2s, "ps2s")
			ps2s.SetOnChange(func() {
				m.gammasdone = false
			})
			ps2s.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			ps2s.SetMin(1e-5)
			ps2s.SetMax(1 - 1e-5)
			ps2s.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(ps2s)
			rs1s := fpg(&m.rs1s, "ln_rs1s")
			rs1s.SetOnChange(func() {
				m.gammasdone = false
			})
			rs1s.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			rs1s.SetMin(-10)
			rs1s.SetMax(-1e-4)
			rs1s.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(rs1s)
			rs2s := fpg(&m.rs2s, "ln_rs2s")
			rs2s.SetOnChange(func() {
				m.gammasdone = false
			})
			rs2s.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			rs2s.SetMin(1e-4)
			rs2s.SetMax(10)
			rs2s.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(rs2s)
		}
	}

	if m.ncatcg > 1 {
		if !m.proportional {
			alphac := fpg(&m.alphac, "alphac")
			alphac.SetOnChange(func() {
				m.gammacdone = false
			})
			alphac.SetPriorFunc(optimize.GammaPrior(1, 2, false))
			alphac.SetMin(1e-2)
			alphac.SetMax(1000)
			alphac.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(alphac)
		} else {
			ps1c := fpg(&m.ps1c, "ps1c")
			ps1c.SetOnChange(func() {
				m.gammacdone = false
			})
			ps1c.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			ps1c.SetMin(1e-5)
			ps1c.SetMax(1 - 1e-5)
			ps1c.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(ps1c)
			ps2c := fpg(&m.ps2c, "ps2c")
			ps2c.SetOnChange(func() {
				m.gammacdone = false
			})
			ps2c.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			ps2c.SetMin(1e-5)
			ps2c.SetMax(1 - 1e-5)
			ps2c.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(ps2c)
			rs1c := fpg(&m.rs1c, "ln_rs1c")
			rs1c.SetOnChange(func() {
				m.gammacdone = false
			})
			rs1c.SetPriorFunc(optimize.UniformPrior(-10, -1e-4, false, false))
			rs1c.SetMin(-10)
			rs1c.SetMax(-1e-4)
			rs1c.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(rs1c)
			rs2c := fpg(&m.rs2c, "ln_rs2c")
			rs2c.SetOnChange(func() {
				m.gammacdone = false
			})
			rs2c.SetPriorFunc(optimize.UniformPrior(1e-4, 10, false, false))
			rs2c.SetMin(1e-4)
			rs2c.SetMax(10)
			rs2c.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(rs2c)
		}
	}
}

// GetParameters returns the model parameter values.
func (m *M8) GetParameters() (p0, p, q, kappa, omega, alphas, alphac float64,
	rs1s, rs2s, ps1s, ps2s float64, rs1c, rs2c, ps1c, ps2c float64) {
	return m.p0, m.p, m.q, m.kappa, m.omega, m.alphas, m.alphac, m.rs1s, m.rs2s, m.ps1s, m.ps2s, m.rs1c, m.rs2c, m.ps1c, m.ps2c
}

// SetParameters sets the model parameter values.
func (m *M8) SetParameters(p0, p, q, kappa, omega, alphas, alphac float64,
	rs1s, rs2s, ps1s, ps2s float64, rs1c, rs2c, ps1c, ps2c float64) {
	if m.addw {
		m.p0 = p0
	} else {
		m.p0 = 1
	}
	m.p = p
	m.q = q
	m.kappa = kappa
	if m.addw && !m.fixw {
		m.omega = omega
	} else {
		m.omega = 1
	}
	m.alphas = alphas
	m.alphac = alphac
	m.gammasdone = false
	m.qbdone = false
	m.q0done = false
}

// SetDefaults sets the default initial parameter values.
func (m *M8) SetDefaults() {
	p0 := 0.69 + rand.Float64()*0.3
	// these initialization values are coming from codeml
	// (p, q, omega); p0 in codeml is 0.9, here I randomize
	// it instead (see above).
	p := 0.2 + rand.Float64()
	q := 1 + rand.Float64()
	omega := 2 + rand.Float64()
	kappa := 1e-2 + rand.Float64()*10
	alphas := 1e-3 + rand.Float64()*10
	alphac := 1e-3 + rand.Float64()*10
	rs1s := -1e-4 - rand.Float64()*(10-1e-4)
	rs2s := 1e-4 + rand.Float64()*(10-1e-4)
	ps1s := 1e-5 + rand.Float64()*(1-2e-5)
	ps2s := 1e-5 + rand.Float64()*(1-2e-5)
	rs1c := -1e-4 - rand.Float64()*(10-1e-4)
	rs2c := 1e-4 + rand.Float64()*(10-1e-4)
	ps1c := 1e-5 + rand.Float64()*(1-2e-5)
	ps2c := 1e-5 + rand.Float64()*(1-2e-5)

	m.SetParameters(p0, p, q, kappa, omega, alphas, alphac,
		rs1s, rs2s, ps1s, ps2s, rs1c, rs2c, ps1c, ps2c)
}

// Organization of the class categories.
// Beta cat 1, internal gamma cat 1, external gamma cat 1
// Beta cat 1, internal gamma cat 1, external gamma cat 2
// ...
// Beta cat 1, internal gamma cat 1, external gamma ncatcg
// ...
// Beta cat 1, internal gamma cat ncatsg^3, external gamma cat 1
// ...
// Beta cat 1, internal gamma cat ncatsg^3, external gamma ncatcg
// Beta cat 2, internal gamma cat 1, external gamma cat 1
// ...
// Beta cat 2, internal gamma cat ncatsg^3, external gamma cat ncatcg
// ...
// ...
// Beta cat ncatb, internal gamma cat 1, external gamma cat 1
// ...
// Beta cat ncatb, internal gamma cat ncatsg^3, external gamma cat ncatcg
// (total: ncatb * ncatsg^3 * ncatcg)
// if m.addw == true [
//   w2, internal gamma cat 1, external gamma cat 1
//   ...
//   w2, internal gamma cat 1, external gamma cat ncatcg
//   w2, internal gamma cat 2, external gamma cat 1
//   ...
//   ...
//   w2, internal gamma cat ncatsg^3, external gamma cat ncatcg
//   (total: (ncatb + 1) * ncatsg^3) * ncatcg
// ]

// UpdateQ updates Q-matrices for after change in the model parameter
// values. This functions updates matrices for the positive selection
// (w2).
func (m *M8) updateQ() {
	gcat := m.ncatsg * m.ncatsg * m.ncatsg

	for c1 := 0; c1 < m.ncatsg; c1++ {
		m.tmp[0] = m.gammas[c1]
		for c2 := 0; c2 < m.ncatsg; c2++ {
			m.tmp[1] = m.gammas[c2]
			for c3 := 0; c3 < m.ncatsg; c3++ {
				m.tmp[2] = m.gammas[c3]

				e := codon.NewEMatrix(m.data.cFreq)

				Q, s := codon.CreateRateTransitionMatrix(m.data.cFreq, m.kappa, m.omega, m.tmp, e.Q)
				e.Set(Q, s)
				err := e.Eigen()
				if err != nil {
					panic("error finding eigen")
				}

				for ecl, rate := range m.gammac {
					catid := (((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl

					e.Copy(m.q0[catid])
					m.q0[catid].ScaleD(rate)

					class := m.ncatb*gcat*m.ncatcg + catid
					for _, node := range m.data.Tree.NodeIDArray() {
						if node == nil {
							continue
						}
						m.qs[class][node.ID] = m.q0[catid]
					}
				}

			}
		}
	}

	m.propdone = false
}

// updateQb updates Q-matrices for after change in the model parameter
// values. This functions updates matrices for the negative selection
// (beta-distributed).
func (m *M8) updateQb() {
	gcat := m.ncatsg * m.ncatsg * m.ncatsg

	// get omega values
	m.omegab = dist.DiscreteBeta(m.p, m.q, m.ncatb, false, m.tmp, m.omegab)

	for c1 := 0; c1 < m.ncatsg; c1++ {
		m.tmp[0] = m.gammas[c1]
		for c2 := 0; c2 < m.ncatsg; c2++ {
			m.tmp[1] = m.gammas[c2]
			for c3 := 0; c3 < m.ncatsg; c3++ {
				m.tmp[2] = m.gammas[c3]

				for icl, omega := range m.omegab {
					e := codon.NewEMatrix(m.data.cFreq)
					Q, s := codon.CreateRateTransitionMatrix(m.data.cFreq, m.kappa, omega, m.tmp, e.Q)
					e.Set(Q, s)
					err := e.Eigen()
					if err != nil {
						panic("error finding eigen")
					}

					for ecl, rate := range m.gammac {
						catid := (icl*gcat+((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl

						e.Copy(m.qb[catid])
						m.qb[catid].ScaleD(rate)

						for _, node := range m.data.Tree.NodeIDArray() {
							if node == nil {
								continue
							}
							m.qs[catid][node.ID] = m.qb[catid]
						}
					}
				}
			}
		}
	}

	m.propdone = false
}

// updateProportions updates proportions if model parameters are
// changing.
func (m *M8) updateProportions() {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	bothcat := m.ncatcg * scat
	p0 := m.p0
	p1 := 1 - p0
	for c1 := 0; c1 < m.ncatsg; c1++ {
		for c2 := 0; c2 < m.ncatsg; c2++ {
			for c3 := 0; c3 < m.ncatsg; c3++ {
				for ecl := range m.gammac {
					catid := (((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl
					pq := m.gammasprop[c1] * m.gammasprop[c2] * m.gammasprop[c3] * m.gammacprop[ecl]
					for i := range m.qb {
						m.prop[0][catid+bothcat*i] = pq * p0 / float64(m.ncatb)
					}
					if m.addw {
						m.prop[0][catid+bothcat*m.ncatb] = pq * p1
					}
				}
			}
		}
	}

	scale := 0.0
	for _, node := range m.data.Tree.NodeIDArray() {
		if node == nil {
			continue
		}
		for i := range m.prop[0] {
			scale += m.prop[0][i] * m.qs[i][node.ID].Scale
		}
		break
	}

	for _, node := range m.data.Tree.NodeIDArray() {
		if node == nil {
			continue
		}
		m.scale[node.ID] = scale
	}

	m.propdone = true
	m.expAllBr = false
}

// cGammaRates returns array with codon gamma rates using NEB approach.
func (m *M8) cGammaRates() []float64 {
	clRates := make([]float64, m.GetNClass())
	scat := m.ncatsg * m.ncatsg * m.ncatsg

	for i := 0; i < m.ncatcg; i++ {
		for j := 0; j < scat; j++ {
			for k := 0; k < m.ncatb; k++ {
				clRates[k*m.ncatcg*scat+j*m.ncatcg+i] = m.gammac[i]
			}
			if m.addw {
				clRates[m.ncatb*m.ncatcg*scat+j*m.ncatcg+i] = m.gammac[i]
			}
		}
	}

	return m.NEBPosterior(clRates)
}

// sGammaRates returns array with site gamma rates using NEB approach.
func (m *M8) sGammaRates() []float64 {
	// these three slices store rates for first, second and third
	// codon positions
	clRates1 := make([]float64, m.GetNClass())
	clRates2 := make([]float64, m.GetNClass())
	clRates3 := make([]float64, m.GetNClass())

	gcat := m.ncatsg * m.ncatsg * m.ncatsg

	for c1 := 0; c1 < m.ncatsg; c1++ {
		rate1 := m.gammas[c1]
		for c2 := 0; c2 < m.ncatsg; c2++ {
			rate2 := m.gammas[c2]
			for c3 := 0; c3 < m.ncatsg; c3++ {
				rate3 := m.gammas[c3]
				for ecl := 0; ecl < m.ncatcg; ecl++ {
					for k := 0; k < m.ncatb; k++ {
						for icl := range m.omegab {
							catid := (icl*gcat+((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl
							clRates1[catid] = rate1
							clRates2[catid] = rate2
							clRates3[catid] = rate3
						}
						if m.addw {
							catid := m.ncatb*gcat*m.ncatcg + (((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl
							clRates1[catid] = rate1
							clRates2[catid] = rate2
							clRates3[catid] = rate3
						}
					}
				}
			}
		}
	}

	// estimate rates for positions separately
	ratePosterior1 := m.NEBPosterior(clRates1)
	ratePosterior2 := m.NEBPosterior(clRates2)
	ratePosterior3 := m.NEBPosterior(clRates3)

	// put rates in the array
	allRatePosterior := make([]float64, m.data.cSeqs.Length()*3)

	for i := range ratePosterior1 {
		allRatePosterior[i*3] = ratePosterior1[i]
		allRatePosterior[i*3+1] = ratePosterior2[i]
		allRatePosterior[i*3+2] = ratePosterior3[i]
	}

	return allRatePosterior
}

// omegaPosterior returns array with omega posterior using NEB approach.
func (m *M8) omegaPosterior() []float64 {
	clOmega := make([]float64, m.GetNClass())
	scat := m.ncatsg * m.ncatsg * m.ncatsg

	for i := 0; i < m.ncatcg; i++ {
		for j := 0; j < scat; j++ {
			for k := 0; k < m.ncatb; k++ {
				clOmega[k*m.ncatcg*scat+j*m.ncatcg+i] = m.omegab[k]
			}
			if m.addw {
				clOmega[m.ncatb*m.ncatcg*scat+j*m.ncatcg+i] = m.omega
			}
		}
	}

	return m.NEBPosterior(clOmega)
}

// Final prints NEB results (only if with positive selection).
func (m *M8) Final(neb, beb, codonRates, siteRates, codonOmega bool) {
	startTime := time.Now()
	defer func() { m.summary.PosteriorTime = time.Since(startTime).Seconds() }()

	if m.ncatcg > 1 && codonRates {
		m.summary.CodonGammaRates = m.cGammaRates()
		log.Notice("Codon gamma rates posterior:", m.summary.CodonGammaRates)
	}

	if m.ncatsg > 1 && siteRates {
		m.summary.SiteGammaRates = m.sGammaRates()
		log.Notice("Site gamma rates posterior:", m.summary.SiteGammaRates)
	}

	if m.ncatb > 1 && codonOmega {
		m.summary.CodonOmega = m.omegaPosterior()
		log.Notice("Codon omega posterior:", m.summary.CodonOmega)
	}

	if neb && m.addw && !m.fixw {
		classes := make([]float64, m.GetNClass())

		gcat := m.ncatsg * m.ncatsg * m.ncatsg

		for c1 := 0; c1 < m.ncatsg; c1++ {
			for c2 := 0; c2 < m.ncatsg; c2++ {
				for c3 := 0; c3 < m.ncatsg; c3++ {
					for ecl := range m.gammac {
						catid := (((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl

						class := m.ncatb*gcat*m.ncatcg + catid
						classes[class] = 1
					}
				}
			}
		}

		posterior := m.NEBPosterior(classes)
		m.summary.SitePosteriorNEB = posterior

		m.PrintPosterior(posterior)
	}
	if m.ncatcg > 1 {
		m.update()
		log.Infof("Codon rates: %v", m.gammac)
		log.Infof("Codon props: %v", m.gammacprop)
	}
	if m.ncatsg > 1 {
		m.update()
		log.Infof("Site rates: %v", m.gammas)
		log.Infof("Site props: %v", m.gammasprop)
	}
}

// update updates matrices and proportions.
func (m *M8) update() {
	if !m.gammasdone {
		if m.ncatsg > 1 {
			if !m.proportional {
				m.gammas = dist.DiscreteGamma(m.alphas, m.alphas, m.ncatsg, false, m.tmp, m.gammas)
			} else {
				rs1s := math.Exp(m.rs1s)
				rs2s := math.Exp(m.rs2s)

				synScale := rs1s*m.ps1s + 1*(1-m.ps1s)*m.ps2s + rs2s*(1-m.ps1s)*(1-m.ps2s)
				m.gammas[0] = rs1s / synScale
				m.gammas[1] = 1 / synScale
				m.gammas[2] = rs2s / synScale
				m.gammasprop[0] = m.ps1s
				m.gammasprop[1] = (1 - m.ps1s) * m.ps2s
				m.gammasprop[2] = (1 - m.ps1s) * (1 - m.ps2s)
			}
			m.q0done = false
			m.qbdone = false
		} else {
			m.gammas[0] = 1
			m.gammasprop[0] = 1
		}
		m.gammasdone = true
	}
	if !m.gammacdone {
		if m.ncatcg > 1 {
			if !m.proportional {
				m.gammac = dist.DiscreteGamma(m.alphac, m.alphac, m.ncatcg, false, m.tmp, m.gammac)
			} else {
				rs1c := math.Exp(m.rs1c)
				rs2c := math.Exp(m.rs2c)

				synScale := rs1c*m.ps1c + 1*(1-m.ps1c)*m.ps2c + rs2c*(1-m.ps1c)*(1-m.ps2c)
				m.gammac[0] = rs1c / synScale
				m.gammac[1] = 1 / synScale
				m.gammac[2] = rs2c / synScale
				m.gammacprop[0] = m.ps1c
				m.gammacprop[1] = (1 - m.ps1c) * m.ps2c
				m.gammacprop[2] = (1 - m.ps1c) * (1 - m.ps2c)
			}
			m.q0done = false
			m.qbdone = false
		} else {
			m.gammac[0] = 1
			m.gammacprop[0] = 1
		}
		m.gammacdone = true
	}
	if !m.qbdone {
		m.updateQb()
	}
	if m.addw && !m.q0done {
		m.updateQ()
	}
	if !m.propdone {
		m.updateProportions()
	}
}

// Summary returns the run summary (site posterior for NEB and BEB).
func (m *M8) Summary() interface{} {
	if !m.summary.Empty() {
		return m.summary
	}
	// nil prevents json from printing "{}"
	return nil
}

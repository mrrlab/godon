package cmodel

import (
	"math"
	"math/rand"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/dist"
	"bitbucket.org/Davydov/godon/optimize"
)

// M0G is an implementation of M0G model.
type M0G struct {
	*BaseModel
	q            []*codon.EMatrix
	omega, kappa float64

	// site gamma alpha parameter
	alphas     float64
	gammas     []float64
	gammasprop []float64
	// parametrization inspired by Scheffler 2006
	// proportions
	ps1s, ps2s float64
	// rates, both are log-scale
	rs1s, rs2s float64

	// codon gamma alpha parameter
	alphac     float64
	gammac     []float64
	gammacprop []float64
	// parametrization inspired by Scheffler 2006
	// proportions
	ps1c, ps2c float64
	// rates, both are log-scale
	rs1c, rs2c float64

	proportional bool

	ncatsg int
	ncatcg int

	qdone      bool
	gammasdone bool
	gammacdone bool

	// temporary array
	tmp []float64
}

// NewM0G creates a new M0G model.
func NewM0G(data *Data, ncatsg, ncatcg int, proportional bool) (m *M0G) {
	gcat := ncatsg * ncatsg * ncatsg

	if proportional {
		if (ncatsg != 3 && ncatsg != 1) || (ncatcg != 3 && ncatcg != 1) {
			log.Fatal("Scheffler's parametrization is only supported for three discrete rates")
		}
	}

	m = &M0G{
		proportional: proportional,
		q:            make([]*codon.EMatrix, gcat*ncatcg),
		gammas:       make([]float64, ncatsg),
		gammasprop:   make([]float64, ncatsg),
		gammac:       make([]float64, ncatcg),
		gammacprop:   make([]float64, ncatcg),
		ncatsg:       ncatsg,
		ncatcg:       ncatcg,
		tmp:          make([]float64, maxInt(ncatsg, ncatcg, 3)),
	}
	for i := 0; i < gcat*ncatcg; i++ {
		m.q[i] = codon.NewEMatrix(data.cFreq)
	}
	if !m.proportional {
		for i := range m.gammacprop {
			m.gammacprop[i] = 1 / float64(len(m.gammacprop))
		}
		for i := range m.gammasprop {
			m.gammasprop[i] = 1 / float64(len(m.gammasprop))
		}
	}

	m.BaseModel = NewBaseModel(data, m)

	m.setupParameters()
	m.setBranchMatrices()
	m.SetDefaults()
	return
}

// GetNClass returns number of site classes.
func (m *M0G) GetNClass() int {
	gcat := m.ncatsg * m.ncatsg * m.ncatsg
	return gcat * m.ncatcg
}

// Copy makes a copy of the model preserving the model parameter
// values.
func (m *M0G) Copy() optimize.Optimizable {
	gcat := m.ncatsg * m.ncatsg * m.ncatsg
	newM := &M0G{
		BaseModel:    m.BaseModel.Copy(),
		proportional: m.proportional,
		q:            make([]*codon.EMatrix, gcat*m.ncatcg),
		gammas:       make([]float64, m.ncatsg),
		gammasprop:   make([]float64, m.ncatsg),
		rs1s:         m.rs1s,
		rs2s:         m.rs2s,
		ps1s:         m.ps1s,
		ps2s:         m.ps2s,
		gammac:       make([]float64, m.ncatcg),
		gammacprop:   make([]float64, m.ncatcg),
		rs1c:         m.rs1c,
		rs2c:         m.rs2c,
		ps1c:         m.ps1c,
		ps2c:         m.ps2c,
		omega:        m.omega,
		kappa:        m.kappa,
		ncatsg:       m.ncatsg,
		ncatcg:       m.ncatcg,
		tmp:          make([]float64, maxInt(m.ncatsg, m.ncatcg, 3)),
	}
	for i := 0; i < gcat*m.ncatcg; i++ {
		newM.q[i] = codon.NewEMatrix(m.data.cFreq)
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
	newM.setBranchMatrices()
	return newM
}

// addParameters adds all the model parameters to the parameter
// storage.
func (m *M0G) addParameters(fpg optimize.FloatParameterGenerator) {
	omega := fpg(&m.omega, "omega")
	omega.SetOnChange(func() {
		m.qdone = false
		m.expAllBr = false
	})
	omega.SetPriorFunc(optimize.GammaPrior(1, 2, false))
	omega.SetProposalFunc(optimize.NormalProposal(0.01))
	omega.SetMin(1e-4)
	omega.SetMax(1000)

	m.parameters.Append(omega)

	kappa := fpg(&m.kappa, "kappa")
	kappa.SetOnChange(func() {
		m.qdone = false
		m.expAllBr = false
	})
	kappa.SetPriorFunc(optimize.UniformPrior(0, 20, false, true))
	kappa.SetProposalFunc(optimize.NormalProposal(0.01))
	kappa.SetMin(1e-2)
	kappa.SetMax(100)

	m.parameters.Append(kappa)

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
			rs1s := fpg(&m.rs1s, "log_rs1s")
			rs1s.SetOnChange(func() {
				m.gammasdone = false
			})
			rs1s.SetPriorFunc(optimize.UniformPrior(-10, -1e-4, false, false))
			rs1s.SetMin(-10)
			rs1s.SetMax(-1e-4)
			rs1s.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(rs1s)
			rs2s := fpg(&m.rs2s, "log_rs2s")
			rs2s.SetOnChange(func() {
				m.gammasdone = false
			})
			rs2s.SetPriorFunc(optimize.UniformPrior(1e-4, 10, false, false))
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
			rs1c := fpg(&m.rs1c, "log_rs1c")
			rs1c.SetOnChange(func() {
				m.gammacdone = false
			})
			rs1c.SetPriorFunc(optimize.UniformPrior(-10, -1e-4, false, false))
			rs1c.SetMin(-10)
			rs1c.SetMax(-1e-4)
			rs1c.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(rs1c)
			rs2c := fpg(&m.rs2c, "log_rs2c")
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
func (m *M0G) GetParameters() (kappa, omega, alphas, alphac float64,
	rs1s, rs2s, ps1s, ps2s float64, rs1c, rs2c, ps1c, ps2c float64) {
	return m.kappa, m.omega, m.alphas, m.alphac, rs1s, rs2s, ps1s, ps2s, rs1c, rs2c, ps1c, ps2c
}

// SetParameters sets the model parameter values.
func (m *M0G) SetParameters(kappa, omega, alphas, alphac float64, rs1s, rs2s, ps1s, ps2s float64, rs1c, rs2c, ps1c, ps2c float64) {
	m.kappa = kappa
	m.omega = omega
	m.alphas = alphas
	m.alphac = alphac
	m.rs1s = rs1s
	m.rs2s = rs2s
	m.ps1s = ps1s
	m.ps2s = ps2s
	m.rs1c = rs1c
	m.rs2c = rs2c
	m.ps1c = ps1c
	m.ps2c = ps2c

	m.qdone = false
	m.gammasdone = false
	m.gammacdone = false
}

// SetDefaults sets the default initial parameter values.
func (m *M0G) SetDefaults() {
	kappa := 1e-2 + rand.Float64()*10
	omega := 0.2 + 0.1*rand.Float64()
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

	m.SetParameters(kappa, omega, alphas, alphac,
		rs1s, rs2s, ps1s, ps2s, rs1c, rs2c, ps1c, ps2c)
}

// Organization of the class categories.
// omega, internal gamma cat 1, external gamma cat 1
// omega, internal gamma cat 1, external gamma cat 2
// ...
// omega, internal gamma cat 1, external gamma ncatcg
// ...
// omega, internal gamma cat ncatsg^3, external gamma cat 1
// ...
// omega, internal gamma cat ncatsg^3, external gamma ncatcg

// setBranchMatrices set matrices for all the branches.
func (m *M0G) setBranchMatrices() {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	for _, node := range m.data.Tree.NodeIDArray() {
		if node == nil {
			continue
		}
		for i := 0; i < m.ncatcg; i++ {
			for j := 0; j < scat; j++ {
				m.qs[i+j*m.ncatcg][node.ID] = m.q[i+j*m.ncatcg]
			}
		}
	}
}

// updateMatrices updates Q-matrices after change in the model parameter
// values.
func (m *M0G) updateMatrices() {
	scale := 0.0
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
					log.Fatal(err)
				}

				for ecl, rate := range m.gammac {
					pq := m.gammasprop[c1] * m.gammasprop[c2] * m.gammasprop[c3] * m.gammacprop[ecl]
					catid := (((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl

					e.Copy(m.q[catid])
					m.q[catid].ScaleD(rate)
					m.prop[0][catid] = pq
					scale += pq * s * rate
				}
			}
		}

	}
	for _, node := range m.data.Tree.NodeIDArray() {
		if node == nil {
			continue
		}
		m.scale[node.ID] = scale
	}

	m.expAllBr = false
}

// update updates matrices and proportions.
func (m *M0G) update() {
	if !m.gammasdone {
		if m.ncatsg > 1 {
			if !m.proportional {
				m.gammas = dist.DiscreteGamma(m.alphas, m.alphas, m.ncatsg, false, m.tmp, m.gammas)
			} else {
				rs1s := math.Exp(m.rs1s)
				rs2s := math.Exp(m.rs2s)

				synScale := rs1s*m.ps1s + 1*(1-m.ps1s)*m.ps2s + 1/rs2s*(1-m.ps1s)*(1-m.ps2s)
				m.gammas[0] = rs1s / synScale
				m.gammas[1] = 1 / synScale
				m.gammas[2] = rs2s / synScale
				m.gammasprop[0] = m.ps1s
				m.gammasprop[1] = (1 - m.ps1s) * m.ps2s
				m.gammasprop[2] = (1 - m.ps1s) * (1 - m.ps2s)
			}
			m.qdone = false
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
			m.qdone = false
		} else {
			m.gammac[0] = 1
			m.gammacprop[0] = 1
		}
		m.gammacdone = true
	}
	if !m.qdone {
		m.updateMatrices()
	}
}

// Final prints codon and site rates for proportional mode.
func (m *M0G) Final(neb, beb, codonRates, siteRates, codonOmega bool) {
	if m.ncatcg > 1 && m.proportional {
		m.update()
		log.Infof("Codon rates: %v", m.gammac)
		log.Infof("Codon props: %v", m.gammacprop)
	}
	if m.ncatsg > 1 && m.proportional {
		m.update()
		log.Infof("Site rates: %v", m.gammas)
		log.Infof("Site props: %v", m.gammasprop)
	}
}

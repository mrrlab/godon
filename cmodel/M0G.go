package cmodel

import (
	"math/rand"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/paml"
)

// M0G is an implementation of M0G model.
type M0G struct {
	*BaseModel
	q            []*codon.EMatrix
	omega, kappa float64

	// site gamma alpha parameter
	alphas float64
	gammas []float64

	// codon gamma alpha parameter
	alphac float64
	gammac []float64

	ncatsg int
	ncatcg int

	qdone      bool
	gammasdone bool
	gammacdone bool

	// temporary array
	tmp []float64
}

// NewM0G creates a new M0G model.
func NewM0G(data *Data, ncatsg, ncatcg int) (m *M0G) {
	gcat := ncatsg * ncatsg * ncatsg
	m = &M0G{
		q:      make([]*codon.EMatrix, gcat*ncatcg),
		gammas: make([]float64, ncatsg),
		gammac: make([]float64, ncatcg),
		ncatsg: ncatsg,
		ncatcg: ncatcg,
		tmp:    make([]float64, maxInt(ncatsg, ncatcg, 3)),
	}
	for i := 0; i < gcat*ncatcg; i++ {
		m.q[i] = codon.NewEMatrix(data.cFreq)
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
		BaseModel: m.BaseModel.Copy(),
		q:         make([]*codon.EMatrix, gcat*m.ncatcg),
		gammas:    make([]float64, m.ncatsg),
		gammac:    make([]float64, m.ncatcg),
		omega:     m.omega,
		kappa:     m.kappa,
		ncatsg:    m.ncatsg,
		ncatcg:    m.ncatcg,
		tmp:       make([]float64, maxInt(m.ncatsg, m.ncatcg, 3)),
	}
	for i := 0; i < gcat*m.ncatcg; i++ {
		newM.q[i] = codon.NewEMatrix(m.data.cFreq)
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
	kappa.SetMax(20)

	m.parameters.Append(kappa)

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
func (m *M0G) GetParameters() (kappa, omega, alphas, alphac float64) {
	return m.kappa, m.omega, m.alphas, m.alphac
}

// SetParameters sets the model parameter values.
func (m *M0G) SetParameters(kappa, omega, alphas, alphac float64) {
	m.kappa = kappa
	m.omega = omega
	m.alphas = alphas
	m.alphac = alphac

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

	m.SetParameters(kappa, omega, alphas, alphac)
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
	gcat := m.ncatsg * m.ncatsg * m.ncatsg
	pq := 1.0 / float64(gcat*m.ncatcg)
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
			m.gammas = paml.DiscreteGamma(m.alphas, m.alphas, m.ncatsg, false, m.tmp, m.gammas)
			m.qdone = false
		} else {
			m.gammas[0] = 1
		}
		m.gammasdone = true
	}
	if !m.gammacdone {
		if m.ncatcg > 1 {
			m.gammac = paml.DiscreteGamma(m.alphac, m.alphac, m.ncatcg, false, m.tmp, m.gammac)
			m.qdone = false
		} else {
			m.gammac[0] = 1
		}
		m.gammacdone = true
	}
	if !m.qdone {
		m.updateMatrices()
	}
}

package cmodel

import (
	"math/rand"

	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/paml"
	"bitbucket.org/Davydov/godon/tree"
)

type M8 struct {
	*BaseModel
	qb           []*EMatrix
	q0           []*EMatrix
	p0           float64
	p, q         float64
	omega, kappa float64
	// site gamma alpha parameter
	alphas float64
	gammas []float64
	// codon gamma alpha parameter
	alphac float64
	gammac []float64
	// omega values for beta distribution
	omegab []float64
	// temporary array for beta computations
	tmp []float64
	// only allow extra omega if addw is true
	addw bool
	// fix omega=1
	fixw       bool
	ncatb      int
	ncatsg     int
	ncatcg     int
	q0done     bool
	qbdone     bool
	propdone   bool
	gammasdone bool
	gammacdone bool
}

func NewM8(cali CodonSequences, t *tree.Tree, cf CodonFrequency, addw, fixw bool, ncatb, ncatsg, ncatcg int) (m *M8) {
	// n site gamma categories, ncatb * n^3 matrices
	gcat := ncatsg * ncatsg * ncatsg
	if ncatb < 2 {
		panic("M8 requires at least two categories")
	}
	m = &M8{
		addw:   addw,
		fixw:   fixw,
		ncatb:  ncatb,
		ncatsg: ncatsg,
		ncatcg: ncatcg,
		qb:     make([]*EMatrix, ncatb*gcat*ncatcg),
		q0:     make([]*EMatrix, gcat*ncatcg),
		gammas: make([]float64, ncatsg),
		gammac: make([]float64, ncatcg),
		tmp:    make([]float64, maxInt(ncatb, ncatsg, ncatcg, 3)),
	}

	for i := 0; i < gcat*ncatcg; i++ {
		m.q0[i] = &EMatrix{}
	}
	for i := 0; i < gcat*ncatcg*ncatb; i++ {
		m.qb[i] = &EMatrix{}
	}

	m.BaseModel = NewBaseModel(cali, t, cf, m)

	m.setupParameters()
	m.SetDefaults()
	return
}

func (m *M8) GetNClass() int {
	gcat := m.ncatsg * m.ncatsg * m.ncatsg
	if m.addw {
		return (gcat * (m.ncatb + 1)) * m.ncatcg
	}
	return gcat * m.ncatb * m.ncatcg
}

func (m *M8) Copy() optimize.Optimizable {
	// n inner gamma categories, ncatb * n^3 matrices
	gcat := m.ncatsg * m.ncatsg * m.ncatsg
	newM := &M8{
		BaseModel: m.BaseModel.Copy(),
		qb:        make([]*EMatrix, m.ncatb*gcat*m.ncatcg),
		q0:        make([]*EMatrix, gcat*m.ncatcg),
		tmp:       make([]float64, maxInt(m.ncatb, m.ncatsg, m.ncatcg, 3)),
		ncatb:     m.ncatb,
		ncatsg:    m.ncatsg,
		ncatcg:    m.ncatcg,
		gammas:    make([]float64, m.ncatsg),
		gammac:    make([]float64, m.ncatcg),
		addw:      m.addw,
		fixw:      m.fixw,
		p0:        m.p0,
		p:         m.p,
		q:         m.q,
		omega:     m.omega,
		kappa:     m.kappa,
		alphas:    m.alphas,
		alphac:    m.alphac,
	}

	for i := 0; i < gcat*m.ncatcg; i++ {
		newM.q0[i] = &EMatrix{}
	}
	for i := 0; i < gcat*m.ncatb*m.ncatcg; i++ {
		newM.qb[i] = &EMatrix{}
	}

	newM.BaseModel.Model = newM
	newM.setupParameters()
	return newM
}

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
	p.SetMin(0.1)
	p.SetMax(100)
	p.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(p)

	q := fpg(&m.q, "q")
	q.SetOnChange(func() {
		m.qbdone = false
	})
	q.SetPriorFunc(optimize.ExponentialPrior(1, false))
	q.SetMin(0.1)
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
	kappa.SetMax(20)
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

func (m *M8) GetParameters() (p0, p, q, kappa, omega, alphas, alphac float64) {
	return m.p0, m.p, m.q, m.kappa, m.omega, m.alphas, m.alphac
}

func (m *M8) SetParameters(p0, p, q, kappa, omega, alphas, alphac float64) {
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

func (m *M8) SetDefaults() {
	p0 := 0.69 + rand.Float64()*0.3
	p := 1e-3 + rand.Float64()*10
	q := 1e-3 + rand.Float64()*50
	kappa := 1e-2 + rand.Float64()*10
	omega := 1.001 + rand.Float64()*10
	alphas := 1e-3 + rand.Float64()*10
	alphac := 1e-3 + rand.Float64()*10
	m.SetParameters(p0, p, q, kappa, omega, alphas, alphac)
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

func (m *M8) updateQ() {
	gcat := m.ncatsg * m.ncatsg * m.ncatsg

	for c1 := 0; c1 < m.ncatsg; c1++ {
		m.tmp[0] = m.gammas[c1]
		for c2 := 0; c2 < m.ncatsg; c2++ {
			m.tmp[1] = m.gammas[c2]
			for c3 := 0; c3 < m.ncatsg; c3++ {
				m.tmp[2] = m.gammas[c3]

				e := &EMatrix{}

				Q, s := createRateTransitionMatrix(m.cf, m.kappa, m.omega, m.tmp, e.Q)
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
					for _, node := range m.tree.NodeIdArray() {
						if node == nil {
							continue
						}
						m.qs[class][node.Id] = m.q0[catid]
					}
				}

			}
		}
	}

	m.propdone = false
}

func (m *M8) updateQb() {
	gcat := m.ncatsg * m.ncatsg * m.ncatsg

	// get omega values
	m.omegab = paml.DiscreteBeta(m.p, m.q, m.ncatb, false, m.tmp, m.omegab)

	for c1 := 0; c1 < m.ncatsg; c1++ {
		m.tmp[0] = m.gammas[c1]
		for c2 := 0; c2 < m.ncatsg; c2++ {
			m.tmp[1] = m.gammas[c2]
			for c3 := 0; c3 < m.ncatsg; c3++ {
				m.tmp[2] = m.gammas[c3]

				for icl, omega := range m.omegab {
					e := &EMatrix{}
					Q, s := createRateTransitionMatrix(m.cf, m.kappa, omega, m.tmp, e.Q)
					e.Set(Q, s)
					err := e.Eigen()
					if err != nil {
						panic("error finding eigen")
					}

					for ecl, rate := range m.gammac {
						catid := (icl*gcat+((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl

						e.Copy(m.qb[catid])
						m.qb[catid].ScaleD(rate)

						for _, node := range m.tree.NodeIdArray() {
							if node == nil {
								continue
							}
							m.qs[catid][node.Id] = m.qb[catid]
						}
					}
				}
			}
		}
	}

	m.propdone = false
}

func (m *M8) updateProportions() {
	scale := 0.0

	gcat := m.ncatsg * m.ncatsg * m.ncatsg

	pqi := m.p0 / float64(m.ncatb*gcat*m.ncatcg)

	if len(m.qb) != m.ncatb*gcat*m.ncatcg {
		panic("wrong number of categories")
	}

	i := 0
	var q *EMatrix
	for i, q = range m.qb {
		m.prop[i] = pqi
		scale += q.Scale * pqi
	}

	if m.addw {
		i++
		pqi := (1 - m.p0) / float64(gcat*m.ncatcg)
		for _, q = range m.q0 {
			scale += q.Scale * pqi
			m.prop[i] = pqi
			i++
		}
		if i != gcat*(m.ncatb+1)*m.ncatcg {
			panic("wrong number of catigories (addw=true)")
		}
	}

	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		m.scale[node.Id] = scale
	}

	m.propdone = true
	m.expAllBr = false
}

func (m *M8) Likelihood() float64 {
	if !m.gammasdone {
		if m.ncatsg > 1 {
			m.gammas = paml.DiscreteGamma(m.alphas, m.alphas, m.ncatsg, false, m.tmp, m.gammas)
			m.q0done = false
			m.qbdone = false
		} else {
			m.gammas[0] = 1
		}
		m.gammasdone = true
	}
	if !m.gammacdone {
		if m.ncatcg > 1 {
			m.gammac = paml.DiscreteGamma(m.alphac, m.alphac, m.ncatcg, false, m.tmp, m.gammac)
			m.q0done = false
			m.qbdone = false
		} else {
			m.gammac[0] = 1
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
	l := m.BaseModel.Likelihood()
	log.Debug("Par:", m.parameters, "L=", l)
	return l
}

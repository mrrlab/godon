package cmodel

import (
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
	// incodon gamma alpha parameter
	alphai float64
	gammai []float64
	// per codon gamma alpha parameter
	alphae float64
	gammae []float64
	// omega values for beta distribution
	omegab []float64
	// temporary array for beta computations
	tmp []float64
	// only allow extra omega if addw is true
	addw bool
	// fix omega=1
	fixw       bool
	ncatb      int
	ncatig     int
	ncateg     int
	q0done     bool
	qbdone     bool
	propdone   bool
	gammaidone bool
	gammaedone bool
}

func NewM8(cali CodonSequences, t *tree.Tree, cf CodonFrequency, addw, fixw bool, ncatb, ncatig, ncateg int) (m *M8) {
	// n inner gamma categories, ncatb * n^3 matrices
	gcat := ncatig * ncatig * ncatig
	if ncatb < 2 {
		panic("M8 requires at least two categories")
	}
	m = &M8{
		addw:   addw,
		fixw:   fixw,
		ncatb:  ncatb,
		ncatig: ncatig,
		ncateg: ncateg,
		qb:     make([]*EMatrix, ncatb*gcat*ncateg),
		q0:     make([]*EMatrix, gcat*ncateg),
		gammai: make([]float64, ncatig),
		gammae: make([]float64, ncateg),
		tmp:    make([]float64, maxInt(ncatb, ncatig, ncateg)),
	}

	for i := 0; i < gcat*ncateg; i++ {
		m.q0[i] = &EMatrix{}
	}
	for i := 0; i < gcat*ncateg*ncatb; i++ {
		m.qb[i] = &EMatrix{}
	}

	m.BaseModel = NewBaseModel(cali, t, cf, m)

	m.setupParameters()
	m.SetDefaults()
	return
}

func (m *M8) GetNClass() int {
	gcat := m.ncatig * m.ncatig * m.ncatig
	if m.addw {
		return (gcat * (m.ncatb + 1)) * m.ncateg
	}
	return gcat * m.ncatb * m.ncateg
}

func (m *M8) Copy() optimize.Optimizable {
	// n inner gamma categories, ncatb * n^3 matrices
	gcat := m.ncatig * m.ncatig * m.ncatig
	newM := &M8{
		BaseModel: m.BaseModel.Copy(),
		qb:        make([]*EMatrix, m.ncatb*gcat*m.ncateg),
		q0:        make([]*EMatrix, gcat*m.ncateg),
		tmp:       make([]float64, maxInt(m.ncatb, m.ncatig, m.ncateg)),
		ncatb:     m.ncatb,
		ncatig:    m.ncatig,
		ncateg:    m.ncateg,
		gammai:    make([]float64, m.ncatig),
		gammae:    make([]float64, m.ncateg),
		addw:      m.addw,
		fixw:      m.fixw,
		p0:        m.p0,
		p:         m.p,
		q:         m.q,
		omega:     m.omega,
		kappa:     m.kappa,
		alphai:    m.alphai,
		alphae:    m.alphae,
	}

	for i := 0; i < gcat*m.ncateg; i++ {
		newM.q0[i] = &EMatrix{}
	}
	for i := 0; i < gcat*m.ncatb*m.ncateg; i++ {
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
	kappa.SetMin(0)
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

	if m.ncatig > 1 {
		alphai := fpg(&m.alphai, "alphai")
		alphai.SetOnChange(func() {
			m.gammaidone = false
		})
		alphai.SetPriorFunc(optimize.GammaPrior(1, 2, false))
		alphai.SetMin(0)
		alphai.SetMax(1000)
		alphai.SetProposalFunc(optimize.NormalProposal(0.01))
		m.parameters.Append(alphai)
	}

	if m.ncateg > 1 {
		alphae := fpg(&m.alphae, "alphae")
		alphae.SetOnChange(func() {
			m.gammaedone = false
		})
		alphae.SetPriorFunc(optimize.GammaPrior(1, 2, false))
		alphae.SetMin(0)
		alphae.SetMax(1000)
		alphae.SetProposalFunc(optimize.NormalProposal(0.01))
		m.parameters.Append(alphae)
	}
}

func (m *M8) GetParameters() (p0, p, q, kappa, omega, alphai, alphae float64) {
	return m.p0, m.p, m.q, m.kappa, m.omega, m.alphai, m.alphae
}

func (m *M8) SetParameters(p0, p, q, kappa, omega, alphai, alphae float64) {
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
	m.alphai = alphai
	m.alphae = alphae
	m.gammaidone = false
	m.qbdone = false
	m.q0done = false
}

func (m *M8) SetDefaults() {
	m.SetParameters(0.5, 2, 2, 2, 0.2, 1, 1)
}

// Organization of the class categories.
// Beta cat 1, internal gamma cat 1, external gamma cat 1
// Beta cat 1, internal gamma cat 1, external gamma cat 2
// ...
// Beta cat 1, internal gamma cat 1, external gamma ncateg
// ...
// Beta cat 1, internal gamma cat ncatig^3, external gamma cat 1
// ...
// Beta cat 1, internal gamma cat ncatig^3, external gamma ncateg
// Beta cat 2, internal gamma cat 1, external gamma cat 1
// ...
// Beta cat 2, internal gamma cat ncatig^3, external gamma cat ncateg
// ...
// ...
// Beta cat ncatb, internal gamma cat 1, external gamma cat 1
// ...
// Beta cat ncatb, internal gamma cat ncatig^3, external gamma cat ncateg
// (total: ncatb * ncatig^3 * ncateg)
// if m.addw == true [
//   w2, internal gamma cat 1, external gamma cat 1
//   ...
//   w2, internal gamma cat 1, external gamma cat ncateg
//   w2, internal gamma cat 2, external gamma cat 1
//   ...
//   ...
//   w2, internal gamma cat ncatig^3, external gamma cat ncateg
//   (total: (ncatb + 1) * ncatig^3) * ncateg
// ]

func (m *M8) updateQ() {
	gcat := m.ncatig * m.ncatig * m.ncatig

	for c1 := 0; c1 < m.ncatig; c1++ {
		m.tmp[0] = m.gammai[c1]
		for c2 := 0; c2 < m.ncatig; c2++ {
			m.tmp[1] = m.gammai[c2]
			for c3 := 0; c3 < m.ncatig; c3++ {
				m.tmp[2] = m.gammai[c3]

				e := &EMatrix{}

				Q, s := createRateTransitionMatrix(m.cf, m.kappa, m.omega, m.tmp, e.Q)
				e.Set(Q, s)
				err := e.Eigen()
				if err != nil {
					panic("error finding eigen")
				}

				for ecl, rate := range m.gammae {
					catid := (((c1*m.ncatig)+c2)*m.ncatig+c3)*m.ncateg + ecl

					e.Copy(m.q0[catid])
					m.q0[catid].ScaleD(rate)

					class := m.ncatb*gcat*m.ncatig + catid
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
	gcat := m.ncatig * m.ncatig * m.ncatig

	// get omega values
	m.omegab = paml.DiscreteBeta(m.p, m.q, m.ncatb, false, m.tmp, m.omegab)

	for c1 := 0; c1 < m.ncatig; c1++ {
		m.tmp[0] = m.gammai[c1]
		for c2 := 0; c2 < m.ncatig; c2++ {
			m.tmp[1] = m.gammai[c2]
			for c3 := 0; c3 < m.ncatig; c3++ {
				m.tmp[2] = m.gammai[c3]

				for icl, omega := range m.omegab {
					e := &EMatrix{}
					Q, s := createRateTransitionMatrix(m.cf, m.kappa, omega, m.tmp, e.Q)
					e.Set(Q, s)
					err := e.Eigen()
					if err != nil {
						panic("error finding eigen")
					}

					for ecl, rate := range m.gammae {
						catid := (icl*gcat+((c1*m.ncatig)+c2)*m.ncatig+c3)*m.ncateg + ecl

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

	gcat := m.ncatig * m.ncatig * m.ncatig

	pqi := m.p0 / float64(m.ncatb*gcat*m.ncateg)

	if len(m.qb) != m.ncatb*gcat*m.ncateg {
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
		pqi := (1 - m.p0) / float64(gcat*m.ncateg)
		for _, q = range m.q0 {
			scale += q.Scale * pqi
			m.prop[i] = pqi
			i++
		}
		if i != gcat*(m.ncatb+1)*m.ncateg {
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
	if !m.gammaidone {
		if m.ncatig > 1 {
			m.gammai = paml.DiscreteGamma(m.alphai, m.alphai, m.ncatig, false, m.tmp, m.gammai)
			m.q0done = false
			m.qbdone = false
		} else {
			m.gammai[0] = 1
		}
		m.gammaidone = true
	}
	if !m.gammaedone {
		if m.ncateg > 1 {
			m.gammae = paml.DiscreteGamma(m.alphae, m.alphae, m.ncateg, false, m.tmp, m.gammae)
			m.q0done = false
			m.qbdone = false
		} else {
			m.gammae[0] = 1
		}
		m.gammaedone = true
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

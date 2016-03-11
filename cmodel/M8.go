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
	q0done     bool
	qbdone     bool
	propdone   bool
	gammaidone bool
}

func NewM8(cali CodonSequences, t *tree.Tree, cf CodonFrequency, addw, fixw bool, ncatb, ncatig int) (m *M8) {
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
		qb:     make([]*EMatrix, ncatb*gcat),
		q0:     make([]*EMatrix, gcat),
		gammai: make([]float64, ncatig),
		tmp:    make([]float64, maxInt(ncatb, ncatig)),
	}

	for i := 0; i < gcat; i++ {
		m.q0[i] = &EMatrix{}
	}
	for i := 0; i < gcat*ncatb; i++ {
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
		return gcat * (m.ncatb + 1)
	}
	return gcat * m.ncatb
}

func (m *M8) Copy() optimize.Optimizable {
	// n inner gamma categories, ncatb * n^3 matrices
	gcat := m.ncatig * m.ncatig * m.ncatig
	newM := &M8{
		BaseModel: m.BaseModel.Copy(),
		qb:        make([]*EMatrix, m.ncatb*gcat),
		q0:        make([]*EMatrix, gcat),
		tmp:       make([]float64, maxInt(m.ncatb, m.ncatig)),
		ncatb:     m.ncatb,
		ncatig:    m.ncatig,
		gammai:    make([]float64, gcat),
		addw:      m.addw,
		fixw:      m.fixw,
		p0:        m.p0,
		p:         m.p,
		q:         m.q,
		omega:     m.omega,
		kappa:     m.kappa,
		alphai:    m.alphai,
	}

	for i := 0; i < gcat; i++ {
		newM.q0[i] = &EMatrix{}
	}
	for i := 0; i < gcat*m.ncatb; i++ {
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
			m.q0done = false
			m.qbdone = false
		})
		alphai.SetPriorFunc(optimize.GammaPrior(1, 2, false))
		alphai.SetMin(0)
		alphai.SetMax(1000)
		alphai.SetProposalFunc(optimize.NormalProposal(0.01))
		m.parameters.Append(alphai)
	}
}

func (m *M8) GetParameters() (p0, p, q, kappa, omega, alphai float64) {
	return m.p0, m.p, m.q, m.kappa, m.omega, m.alphai
}

func (m *M8) SetParameters(p0, p, q, kappa, omega, alphai float64) {
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
	m.gammaidone = false
	m.qbdone = false
	m.q0done = false
}

func (m *M8) SetDefaults() {
	m.SetParameters(0.5, 2, 2, 2, 0.2, 1)
}

// Organization of the class categories.
// Beta cat 1, internal gamma cat 1
// Beta cat 1, internal gamma cat 2
// ...
// Beta cat 1, internal gamma cat ncatig^3
// Beta cat 2, internal gamma cat 1
// Beta cat 2, internal gamma cat 2
// ...
// Beta cat 2, internal gamma cat ncatig^3
// ...
// ...
// Beta cat ncatb, internal gamma cat 1
// Beta cat ncatb, internal gamma cat 2
// ...
// Beta cat ncatb, internal gamma cat ncatig^3
// (total: ncatb * ncatig^3)
// if m.addw == true [
//   w2, internal gamma cat 1
//   w2, internal gamma cat 2
//   ...
//   w2, internal gamma cat ncatig^3
//   (total: (ncatb + 1) * ncatig^3)
// ]

func (m *M8) updateQ() {
	gcat := m.ncatig * m.ncatig * m.ncatig

	for c1 := 0; c1 < m.ncatig; c1++ {
		m.tmp[0] = m.gammai[c1]
		for c2 := 0; c2 < m.ncatig; c2++ {
			m.tmp[1] = m.gammai[c2]
			for c3 := 0; c3 < m.ncatig; c3++ {
				m.tmp[2] = m.gammai[c3]

				catid := ((c1*m.ncatig)+c2)*m.ncatig + c3

				Q, s := createRateTransitionMatrix(m.cf, m.kappa, m.omega, m.tmp, m.q0[catid].Q)
				m.q0[catid].Set(Q, s)
				err := m.q0[catid].Eigen()
				if err != nil {
					panic("error finding eigen")
				}

				class := m.ncatb*gcat + catid
				for _, node := range m.tree.NodeIdArray() {
					if node == nil {
						continue
					}
					m.qs[class][node.Id] = m.q0[catid]
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

				catid := ((c1*m.ncatig)+c2)*m.ncatig + c3

				for i, omega := range m.omegab {
					ind := i*gcat + catid
					Q, s := createRateTransitionMatrix(m.cf, m.kappa, omega, m.tmp, m.qb[ind].Q)
					m.qb[ind].Set(Q, s)

					err := m.qb[ind].Eigen()
					if err != nil {
						panic("error finding eigen")
					}
					for _, node := range m.tree.NodeIdArray() {
						if node == nil {
							continue
						}
						m.qs[ind][node.Id] = m.qb[ind]
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

	pqi := m.p0 / float64(m.ncatb*gcat)

	if len(m.qb) != m.ncatb*gcat {
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
		pqi := (1 - m.p0) / float64(gcat)
		for _, q = range m.q0 {
			scale += q.Scale * pqi
			m.prop[i] = pqi
			i++
		}
		if i != gcat*(m.ncatb+1) {
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
			m.gammaidone = true
			m.q0done = false
			m.qbdone = false
		} else {
			m.gammai[0] = 1
		}
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
	return m.BaseModel.Likelihood()
}

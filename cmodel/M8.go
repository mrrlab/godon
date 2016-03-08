package cmodel

import (
	//	"fmt"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/paml"
	"bitbucket.org/Davydov/godon/tree"
)

type M8 struct {
	*BaseModel
	qb           []*EMatrix
	q0           *EMatrix
	p0           float64
	p, q         float64
	omega, kappa float64
	// omega values for beta distribution
	omegab []float64
	// temporary array for beta computations
	tmp      []float64
	ncat     int
	q0done   bool
	qbdone   bool
	propdone bool
}

func NewM8(cali CodonSequences, t *tree.Tree, cf CodonFrequency, ncat int) (m *M8) {
	m = &M8{
		ncat: ncat,
		qb:   make([]*EMatrix, ncat),
		q0:   &EMatrix{},
		tmp:  make([]float64, ncat),
	}

	for i := 0; i < ncat; i++ {
		m.qb[i] = &EMatrix{}
	}

	m.BaseModel = NewBaseModel(cali, t, cf, m)

	m.setupParameters()
	m.SetDefaults()
	return
}

func (m *M8) GetNClass() int {
	return m.ncat + 1
}

func (m *M8) Copy() optimize.Optimizable {
	newM := &M8{
		BaseModel: m.BaseModel.Copy(),
		qb:        make([]*EMatrix, m.ncat),
		tmp:       make([]float64, m.ncat),
		ncat:      m.ncat,
		q0:        &EMatrix{},
		p0:        m.p0,
		p:         m.p,
		q:         m.q,
		omega:     m.omega,
		kappa:     m.kappa,
	}

	for i := 0; i < m.ncat; i++ {
		newM.qb[i] = &EMatrix{}
	}

	newM.BaseModel.Model = newM
	newM.setupParameters()
	return newM
}

func (m *M8) addParameters(fpg optimize.FloatParameterGenerator) {
	p0 := fpg(&m.p0, "p0")
	p0.SetOnChange(func() {
		m.propdone = false
	})
	p0.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
	p0.SetMin(0)
	p0.SetMax(1)
	p0.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(p0)

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

func (m *M8) GetParameters() (p0, p, q, kappa, omega float64) {
	return m.p0, m.p, m.q, m.kappa, m.omega
}

func (m *M8) SetParameters(p0, p, q, kappa, omega float64) {
	m.p0 = p0
	m.p = p
	m.q = q
	m.kappa = kappa
	m.omega = omega
	m.qbdone = false
	m.q0done = false
}

func (m *M8) SetDefaults() {
	m.SetParameters(0.5, 2, 2, 2, 0.2)
}

func (m *M8) updateQ() {
	Q, s := createTransitionMatrix(m.cf, m.kappa, m.omega, m.q0.Q)
	m.q0.Set(Q, s)

	err := m.q0.Eigen()
	if err != nil {
		panic("error finding eigen")
	}
	class := m.ncat
	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		m.qs[class][node.Id] = m.q0
	}

	m.propdone = false
}

func (m *M8) updateQb() {
	// first get omega values
	m.omegab = paml.DiscreteBeta(m.p, m.q, m.ncat, false, m.tmp, m.omegab)
	//fmt.Println(m.p, m.q, m.omegab)

	for i, omega := range m.omegab {
		Q, s := createTransitionMatrix(m.cf, m.kappa, omega, m.qb[i].Q)
		m.qb[i].Set(Q, s)

		err := m.qb[i].Eigen()
		if err != nil {
			panic("error finding eigen")
		}
		for _, node := range m.tree.NodeIdArray() {
			if node == nil {
				continue
			}
			m.qs[i][node.Id] = m.qb[i]
		}

	}

	m.propdone = false
}

func (m *M8) updateProportions() {
	scale := m.q0.Scale * (1 - m.p0)
	m.prop[m.ncat] = 1 - m.p0
	pqi := m.p0 / float64(m.ncat)
	for i, q := range m.qb {
		m.prop[i] = pqi
		scale += q.Scale * pqi
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
	if !m.qbdone {
		m.updateQb()
	}
	if !m.q0done {
		m.updateQ()
	}
	if !m.propdone {
		m.updateProportions()
	}
	return m.BaseModel.Likelihood()
}

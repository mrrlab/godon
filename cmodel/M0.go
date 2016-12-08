package cmodel

import (
	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
)

// M0 is an implementation of M0 model.
type M0 struct {
	*BaseModel
	q            *codon.EMatrix
	omega, kappa float64
	qdone        bool
}

// NewM0 creates a new M0 model.
func NewM0(data *Data) (m *M0) {
	m = &M0{
		q: codon.NewEMatrix(data.cFreq),
	}
	m.BaseModel = NewBaseModel(data, m)
	m.prop[0][0] = 1

	m.setupParameters()
	m.SetDefaults()
	return
}

// GetNClass returns number of site classes.
func (m *M0) GetNClass() int {
	return 1
}

// Copy makes a copy of the model preserving the model parameter
// values.
func (m *M0) Copy() optimize.Optimizable {
	newM := &M0{
		BaseModel: m.BaseModel.Copy(),
		q:         codon.NewEMatrix(m.data.cFreq),
		omega:     m.omega,
		kappa:     m.kappa,
	}
	newM.BaseModel.model = newM
	newM.setupParameters()
	return newM
}

// addParameters adds all the model parameters to the parameter
// storage.
func (m *M0) addParameters(fpg optimize.FloatParameterGenerator) {
	omega := fpg(&m.omega, "omega")
	omega.SetOnChange(func() {
		m.qdone = false
		m.expAllBr = false
	})
	omega.SetPriorFunc(optimize.GammaPrior(1, 2, false))
	omega.SetProposalFunc(optimize.NormalProposal(0.01))
	omega.SetMin(1e-4)
	omega.SetMax(1000)

	kappa := fpg(&m.kappa, "kappa")
	kappa.SetOnChange(func() {
		m.qdone = false
		m.expAllBr = false
	})
	kappa.SetPriorFunc(optimize.UniformPrior(0, 20, false, true))
	kappa.SetProposalFunc(optimize.NormalProposal(0.01))
	kappa.SetMin(1e-2)
	kappa.SetMax(20)

	m.parameters.Append(omega)
	m.parameters.Append(kappa)
}

// GetParameters returns the model parameter values.
func (m *M0) GetParameters() (kappa, omega float64) {
	return m.kappa, m.omega
}

// SetParameters sets the model parameter values.
func (m *M0) SetParameters(kappa, omega float64) {
	m.kappa = kappa
	m.omega = omega
	m.qdone = false
}

// SetDefaults sets the default initial parameter values.
func (m *M0) SetDefaults() {
	m.SetParameters(1, 1)
}

// UpdateMatrix updates Q-matrix after change in the model parameter
// values.
func (m *M0) UpdateMatrix() {
	Q, s := codon.CreateTransitionMatrix(m.data.cFreq, m.kappa, m.omega, m.q.Q)
	m.q.Set(Q, s)

	err := m.q.Eigen()
	if err != nil {
		panic("error finding eigen")
	}
	for i := 0; i < len(m.qs[0]); i++ {
		m.qs[0][i] = m.q
		m.scale[i] = m.q.Scale
	}
	m.expAllBr = false
}

// update updates matrices and proportions.
func (m *M0) update() {
	if !m.qdone {
		m.UpdateMatrix()
	}
}

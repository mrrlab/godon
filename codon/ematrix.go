package codon

import (
	"errors"
	"math"

	"github.com/skelterjohn/go.matrix"
)

const smallFreq = 1e-20

type EMatrix struct {
	Q     *matrix.DenseMatrix
	Scale float64
	CF    CodonFrequency
	v     *matrix.DenseMatrix
	d     *matrix.DenseMatrix
	iv    *matrix.DenseMatrix
}

func NewEMatrix(Q *matrix.DenseMatrix, scale float64) *EMatrix {
	return &EMatrix{Q: Q, Scale: scale}
}

func (m *EMatrix) Copy() (newM *EMatrix) {
	newM = &EMatrix{Scale: m.Scale}
	if m.Q != nil {
		newM.Q = m.Q
	}
	if m.v != nil {
		newM.v = m.v
	}
	if m.d != nil {
		newM.d = m.d
	}
	if m.iv != nil {
		newM.iv = m.iv
	}
	return
}

func (m *EMatrix) Set(Q *matrix.DenseMatrix, scale float64) {
	m.Q = Q
	m.Scale = scale
	m.v = nil
}

func (m *EMatrix) Eigen() (err error) {
	// Please refer to the EigenQREV pdf from PAML for explanations.
	if m.v != nil {
		return nil
	}
	if m.Scale < smallScale {
		// No need to eigen 0-matrix
		return nil
	}
	rows := m.Q.Rows()
	cols := m.Q.Cols()

	// make sure we have no zero frequences
	psum := 0.0
	if len(m.CF) == 0 {
		return errors.New("EMatrix has an empty codon frequency")
	}
	for _, p := range m.CF {
		p = math.Max(smallFreq, p)
		psum += p
	}

	Pi := matrix.Zeros(cols, rows)
	Pi_i := matrix.Zeros(cols, rows)
	for i, p := range m.CF {
		p = math.Max(smallFreq, p) / psum
		Pi.Set(i, i, math.Sqrt(p))
		Pi_i.Set(i, i, 1/math.Sqrt(p))
	}

	// Compute symmetric matrix A = Pi * Q * Pi_i
	A := matrix.Product(Pi, m.Q, Pi_i)

	R, D, err := A.Eigen()
	if err != nil {
		return err
	}

	m.d = D

	m.v = matrix.Product(Pi_i, R)

	err = R.TransposeInPlace()
	if err != nil {
		return err
	}

	m.iv = matrix.Product(R, Pi)
	return nil
}

func (m *EMatrix) Exp(cD *matrix.DenseMatrix, t float64) (*matrix.DenseMatrix, error) {
	if m.d.Cols() != m.d.Rows() {
		return nil, errors.New("D isn't a square matrix")
	}
	if m.Scale < smallScale {
		return IdentityP, nil
	}
	// This is a dirty hack to allow 0-scale matricies
	if math.IsInf(t, 1) {
		t = math.MaxFloat64
	}

	for i := 0; i < m.d.Rows(); i++ {
		cD.Set(i, i, math.Exp(m.d.Get(i, i)*t))
	}
	return matrix.Product(m.v, cD, m.iv), nil
}

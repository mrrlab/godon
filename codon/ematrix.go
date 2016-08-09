package codon

import (
	"errors"
	"math"

	"github.com/skelterjohn/go.matrix"
)

// EMatrix stores Q-matrix and it's eigendecomposition to quickly
// compute e^Qt.
type EMatrix struct {
	// Q is Q-matrix
	Q     *matrix.DenseMatrix
	// Scale is matrix scale.
	Scale float64
	// CF is codon frequency.
	CF    CodonFrequency
	v     *matrix.DenseMatrix
	d     *matrix.DenseMatrix
	iv    *matrix.DenseMatrix
}

// NewEMatrix creates a new EMatrix.
func NewEMatrix(Q *matrix.DenseMatrix, scale float64) *EMatrix {
	return &EMatrix{Q: Q, Scale: scale}
}

// Copy creates a copy of EMatrix while saving eigendecomposition.
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

// Set sets Q-matrix and its' scale.
func (m *EMatrix) Set(Q *matrix.DenseMatrix, scale float64) {
	m.Q = Q
	m.Scale = scale
	m.v = nil
}

// Eigen performs eigendecomposition.
func (m *EMatrix) Eigen() (err error) {
	if m.v != nil {
		return nil
	}
	m.v, m.d, err = m.Q.Eigen()
	if err != nil {
		return err
	}
	m.iv, err = m.v.Inverse()
	if err != nil {
		return err
	}
	return nil
}

// Exp computes P=e^Qt and writes it to cD matrix.
func (m *EMatrix) Exp(cD *matrix.DenseMatrix, t float64) (*matrix.DenseMatrix, error) {
	if m.d.Cols() != m.d.Rows() {
		return nil, errors.New("D isn't a square matrix")
	}
	if m.Scale*t < smallScale {
		return IdentityP, nil
	}
	// This is a dirty hack to allow 0-scale matricies
	if math.IsInf(t, 1) {
		t = math.MaxFloat64
	}

	for i := 0; i < m.d.Rows(); i++ {
		cD.Set(i, i, math.Exp(m.d.Get(i, i)*t))
	}
	P := matrix.Product(m.v, cD, m.iv)
	Parr := P.Array()
	// Make sure there's no negative elements
	for i := range Parr {
		Parr[i] = math.Abs(Parr[i])
	}
	return P, nil
}

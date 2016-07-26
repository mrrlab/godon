package codon

import (
	"errors"
	"math"

	"github.com/skelterjohn/go.matrix"
)

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

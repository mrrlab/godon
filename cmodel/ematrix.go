package cmodel

import (
	"errors"
	"math"

	"github.com/skelterjohn/go.matrix"
)

type EMatrix struct {
	Q     *matrix.DenseMatrix
	Scale float64
	v     *matrix.DenseMatrix
	d     *matrix.DenseMatrix
	iv    *matrix.DenseMatrix
}

func NewEMatrix(Q *matrix.DenseMatrix, scale float64) *EMatrix {
	return &EMatrix{Q: Q, Scale: scale}
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
	for i := 0; i < m.d.Rows(); i++ {
		cD.Set(i, i, math.Exp(m.d.Get(i, i)*t))
	}
	return matrix.Product(m.v, cD, m.iv), nil
}

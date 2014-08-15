package main

import (
	"errors"
	"math"

	"github.com/skelterjohn/go.matrix"
)

type EMatrix struct {
	Q  *matrix.DenseMatrix
	v  *matrix.DenseMatrix
	d  *matrix.DenseMatrix
	iv *matrix.DenseMatrix
}

func NewEMatrix(Q *matrix.DenseMatrix) *EMatrix {
	return &EMatrix{Q, nil, nil, nil}
}

func (m *EMatrix) Set(Q *matrix.DenseMatrix) {
	m.Q = Q
	m.v = nil
	m.d = nil
	m.iv = nil
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

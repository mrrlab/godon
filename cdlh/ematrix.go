package main

import (
	"errors"
	"math"

	"github.com/skelterjohn/go.matrix"
)

type EMatrix struct {
	Q  *matrix.DenseMatrix
	V  *matrix.DenseMatrix
	D  *matrix.DenseMatrix
	iV *matrix.DenseMatrix
}

func NewEMatrix(Q *matrix.DenseMatrix) *EMatrix {
	return &EMatrix{Q, nil, nil, nil}
}

func (m *EMatrix) Eigen() (err error) {
	if m.V != nil {
		return nil
	}
	m.V, m.D, err = m.Q.Eigen()
	if err != nil {
		return err
	}
	m.iV, err = m.V.Inverse()
	if err != nil {
		return err
	}
	return nil
}

func (m *EMatrix) Exp(cD *matrix.DenseMatrix, t float64) (*matrix.DenseMatrix, error) {
	if m.D.Cols() != m.D.Rows() {
		return nil, errors.New("D isn't a square matrix")
	}
	for i := 0; i < m.D.Rows(); i++ {
		cD.Set(i, i, math.Exp(m.D.Get(i, i)*t))
	}
	return matrix.Product(m.V, cD, m.iV), nil
}

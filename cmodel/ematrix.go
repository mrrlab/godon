package cmodel

import (
	"errors"
	"math"

	"github.com/gonum/matrix/mat64"
)

type EMatrix struct {
	Q     *mat64.Dense
	Scale float64
	v     *mat64.Dense
	d     *mat64.Dense
	iv    *mat64.Dense
}

func NewEMatrix(Q *mat64.Dense, scale float64) *EMatrix {
	//cols, rows := Q.Dims()
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

func (m *EMatrix) Set(Q *mat64.Dense, scale float64) {
	m.Q = Q
	m.Scale = scale
	m.v = nil
}

func (m *EMatrix) ScaleD(scale float64) {
	if m.d == nil {
		panic("Scaling a nil matrix")
	}
	m.d = scaleMatrix(m.d, scale, nil)
	m.Scale *= scale
}

func (m *EMatrix) Eigen() (err error) {
	if m.v != nil {
		return nil
	}
	rows, cols := m.Q.Dims()
	if m.iv == nil {
		m.iv = mat64.NewDense(cols, rows, nil)
	}

	decomp := mat64.Eigen(m.Q, 1E-8)
	m.v = decomp.V
	m.d = decomp.D()
	err = m.iv.Inverse(m.v)
	if err != nil {
		return err
	}
	return nil
}

func (m *EMatrix) Exp(cD *mat64.Dense, t float64) (*mat64.Dense, error) {
	rows, cols := m.Q.Dims()
	if cols != rows {
		return nil, errors.New("D isn't a square matrix")
	}
	// This is a dirty hack to allow 0-scale matricies
	if math.IsInf(t, 1) {
		t = math.MaxFloat64
	}

	for i := 0; i < rows; i++ {
		cD.Set(i, i, math.Exp(m.d.At(i, i)*t))
	}
	res := mat64.NewDense(cols, rows, nil)
	res.Mul(m.v, cD)
	res.Mul(res, m.iv)
	// Remove sligtly negative values
	res.Apply(func(r, c int, v float64) float64 {
		return math.Max(0, v)
	}, res)
	return res, nil
}

package codon

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

func (m *EMatrix) Copy(recv *EMatrix) *EMatrix {
	if recv == nil {
		recv = &EMatrix{}
	}
	recv.Q = m.Q
	recv.v = m.v
	recv.d = m.d
	recv.iv = m.iv
	recv.Scale = m.Scale
	return recv
}

func (m *EMatrix) Set(Q *mat64.Dense, scale float64) {
	m.Q = Q
	m.Scale = scale
	m.v = nil
}

func (m *EMatrix) ScaleD(scale float64) {
	if m.Scale < smallScale {
		// no need to scale almost zero matrix
		return
	}
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
	if m.Scale < smallScale {
		// No need to eigen 0-matrix
		return nil
	}
	rows, cols := m.Q.Dims()
	if m.iv == nil {
		m.iv = mat64.NewDense(cols, rows, nil)
	}

	decomp := mat64.Eigen{}

	status := decomp.Factorize(m.Q, true)

	if !status {
		panic("Error decompozing Q")
	}
	m.v = decomp.Vectors()
	d := decomp.Values(nil)
	m.d = mat64.NewDense(cols, rows, nil)

	for i, v := range d {
		m.d.Set(i, i, real(v))
	}
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
	if m.Scale < smallScale {
		return IdentityP, nil
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

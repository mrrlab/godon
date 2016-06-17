package codon

import (
	"errors"
	"math"

	"github.com/gonum/matrix/mat64"
)

const smallFreq = 1e-20

type EMatrix struct {
	Q     *mat64.Dense
	Scale float64
	CF    CodonFrequency
	v     *mat64.Dense
	d     *mat64.Dense
	iv    *mat64.Dense
}

func NewEMatrix(Q *mat64.Dense, scale float64, cf CodonFrequency) *EMatrix {
	//cols, rows := Q.Dims()
	return &EMatrix{Q: Q, Scale: scale, CF: cf}
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
	recv.CF = m.CF
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
	// Please refer to the EigenQREV pdf from PAML for explanations.
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

	// make sure we have no zero frequences
	psum := 0.0
	for _, p := range m.CF {
		p = math.Max(smallFreq, p)
		psum += p
	}

	// First compute matrix Pi=p_1^{1/2}, p_2^{1/2}. ...
	// and Pi_i (inverse)
	Pi := mat64.NewDense(cols, rows, nil)
	Pi_i := mat64.NewDense(cols, rows, nil)
	for i, p := range m.CF {
		p = math.Max(smallFreq, p) / psum
		Pi.Set(i, i, math.Sqrt(p))
		Pi_i.Set(i, i, 1/math.Sqrt(p))
	}

	// Compute symmetric matrix A = Pi * Q * Pi_i
	A := mat64.NewDense(cols, rows, nil)
	A.Mul(Pi, m.Q)
	A.Mul(A, Pi_i)
	AS := mat64.NewSymDense(cols, A.RawMatrix().Data)

	decomp := mat64.EigenSym{}

	status := decomp.Factorize(AS, true)

	if !status {
		panic("Error decompozing Q")
	}
	R := mat64.NewDense(cols, rows, nil)
	R.EigenvectorsSym(&decomp)
	m.v = mat64.NewDense(cols, rows, nil)
	m.v.Mul(Pi_i, R)

	d := decomp.Values(nil)
	m.d = mat64.NewDense(cols, rows, nil)

	for i, v := range d {
		m.d.Set(i, i, v)
	}
	err = m.iv.Inverse(R)
	if err != nil {
		return err
	}
	m.iv.Mul(m.iv, Pi)
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

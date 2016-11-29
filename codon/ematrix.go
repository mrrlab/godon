package codon

import (
	"errors"
	"math"

	"github.com/gonum/matrix/mat64"
)

// smallFreq is a small frequency which is used instead of zeros in
// the codon frequency.
const smallFreq = 1e-20

// EMatrix stores Q-matrix and it's eigendecomposition to quickly
// compute e^Qt.
type EMatrix struct {
	// Q is Q-matrix
	Q *mat64.Dense
	// Scale is matrix scale.
	Scale float64
	// CF is codon frequency.
	CF Frequency
	v  *mat64.Dense
	d  *mat64.Dense
	iv *mat64.Dense
}

// NewEMatrix creates a new EMatrix.
func NewEMatrix(cf Frequency) *EMatrix {
	//cols, rows := Q.Dims()
	return &EMatrix{CF: cf}
}

// Copy creates a copy of EMatrix while saving eigendecomposition.
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

// Set sets Q-matrix and its' scale.
func (m *EMatrix) Set(Q *mat64.Dense, scale float64) {
	m.Q = Q
	m.Scale = scale
	m.v = nil
}

// ScaleD scales matrix after the eigendecomposition.
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

// Eigen performs eigendecomposition.
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
	if len(m.CF.Freq) == 0 {
		return errors.New("EMatrix has an empty codon frequency")
	}
	for _, p := range m.CF.Freq {
		p = math.Max(smallFreq, p)
		psum += p
	}

	// First compute matrix Pi=p_1^{1/2}, p_2^{1/2}. ...
	// and Pi_i (inverse)
	Pi := mat64.NewDense(cols, rows, nil)
	PiI := mat64.NewDense(cols, rows, nil)
	for i, p := range m.CF.Freq {
		p = math.Max(smallFreq, p) / psum
		Pi.Set(i, i, math.Sqrt(p))
		PiI.Set(i, i, 1/math.Sqrt(p))
	}

	// Compute symmetric matrix A = Pi * Q * Pi_i
	A := mat64.NewDense(cols, rows, nil)
	A.Mul(Pi, m.Q)
	A.Mul(A, PiI)
	AS := mat64.NewSymDense(cols, A.RawMatrix().Data)

	decomp := mat64.EigenSym{}

	status := decomp.Factorize(AS, true)

	if !status {
		panic("Error decompozing Q")
	}
	R := mat64.NewDense(cols, rows, nil)
	R.EigenvectorsSym(&decomp)
	m.v = mat64.NewDense(cols, rows, nil)
	m.v.Mul(PiI, R)

	d := decomp.Values(nil)
	m.d = mat64.NewDense(cols, rows, nil)

	for i, v := range d {
		m.d.Set(i, i, v)
	}

	m.iv.Mul(R.T(), Pi)
	return nil
}

// Exp computes P=e^Qt and writes it to cD matrix.
func (m *EMatrix) Exp(cD *mat64.Dense, t float64) (*mat64.Dense, error) {
	rows, cols := m.Q.Dims()
	if cols != rows {
		return nil, errors.New("D isn't a square matrix")
	}
	if m.Scale*t < smallScale {
		return createIdentityMatrix(cols), nil
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
	rawRes := res.RawMatrix().Data
	for i := range rawRes {
		if rawRes[i] < 0 {
			rawRes[i] = 0
		}
	}
	return res, nil
}

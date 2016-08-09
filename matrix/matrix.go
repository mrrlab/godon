// Matrix is a simplistic wrapper around GSL matrix type.
package matrix

// #cgo pkg-config: gsl
// #include <gsl/gsl_linalg.h>
import "C"
import (
	"bytes"
	"errors"
	"strconv"
)

// Matrix is a matrix.
type Matrix struct {
	data       []float64
	matrixView C.gsl_matrix_view
}

// NewFromArray creates a new matrix of size n1*n2 from an array.
func NewFromArray(data []float64, n1 int, n2 int) (*Matrix, error) {
	if len(data) != n1*n2 {
		return nil, errors.New("matrix dimensions don't match slice size")
	}
	return &Matrix{data, C.gsl_matrix_view_array((*C.double)(&data[0]),
		C.size_t(n1), C.size_t(n2))}, nil

}

// New creates a new matrix.
func New(n1 int, n2 int) (*Matrix, error) {
	if n1 < 1 || n2 < 1 {
		return nil, errors.New("matrix dimensions should be > 0")
	}
	data := make([]float64, n1*n2)
	return &Matrix{data, C.gsl_matrix_view_array((*C.double)(&data[0]),
		C.size_t(n1), C.size_t(n2))}, nil
}

// String returns matrix in printable form.
func (m *Matrix) String() string {
	var buffer bytes.Buffer
	if m.data == nil {
		return "<Uninitialized matrix>"
	}
	buffer.WriteString("<Matrix\n")
	for i1 := C.size_t(0); i1 < m.matrixView.matrix.size1; i1++ {
		if i1 == 10 {
			buffer.WriteString("...\n")
			break
		}
		buffer.WriteString("  ")
		for i2 := C.size_t(0); i2 < m.matrixView.matrix.size2; i2++ {
			if i2 == 10 {
				buffer.WriteString("...")
				break
			}
			buffer.WriteString(strconv.FormatFloat(m.data[int(i1)*int(m.matrixView.matrix.size2)+int(i2)], 'E', 3, 64))
			if i2 < m.matrixView.matrix.size2-1 {
				buffer.WriteByte('\t')
			}
		}
		buffer.WriteByte('\n')
	}
	buffer.WriteByte('>')
	return buffer.String()
}

// Empty creates a new empty matrix with a size of the old one.
func (m *Matrix) Empty() *Matrix {
	nm, _ := New(
		int(m.matrixView.matrix.size1),
		int(m.matrixView.matrix.size2))
	return nm
}

// Copy copies content of matrix into dest.
func (m *Matrix) Copy(dest *Matrix) {
	C.gsl_matrix_memcpy(&dest.matrixView.matrix, &m.matrixView.matrix)
}

// Scale scales a matrix.
func (m *Matrix) Scale(x float64) {
	C.gsl_matrix_scale(&m.matrixView.matrix, C.double(x))
}

// Exponential computes a matrix exponential.
func (m *Matrix) Exponential(em *Matrix) {
	C.gsl_linalg_exponential_ss(&m.matrixView.matrix, &em.matrixView.matrix, 0)
}

// SetItem sets an element of a matrix.
func (m *Matrix) SetItem(i1, i2 int, x float64) {
	m.data[i1*int(m.matrixView.matrix.size2)+i2] = x
}

// ScaleItem scales an element of a matrix.
func (m *Matrix) ScaleItem(i1, i2 int, x float64) {
	m.data[i1*int(m.matrixView.matrix.size2)+i2] = m.data[i1*int(m.matrixView.matrix.size2)+i2] * x
}

// GetItem returns an element of a matrix.
func (m *Matrix) GetItem(i1, i2 int) float64 {
	return m.data[i1*int(m.matrixView.matrix.size2)+i2]
}

// GetSize returns matrix dimensions.
func (m *Matrix) GetSize() (int, int) {
	return int(m.matrixView.matrix.size1), int(m.matrixView.matrix.size2)
}

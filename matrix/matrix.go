package matrix

// #cgo pkg-config: gsl
// #include <gsl/gsl_linalg.h>
import "C"
import (
	"errors"
	"strconv"
	"bytes"
)

type Matrix struct {
	data []float64
	matrixView C.gsl_matrix_view
}

func NewFromArray(data []float64, n1 int, n2 int) (*Matrix, error) {
	if len(data) != n1 * n2  {
		return nil, errors.New("matrix dimensions don't match slice size")
	}
	return &Matrix{data, C.gsl_matrix_view_array((*C.double)(&data[0]),
		C.size_t(n1), C.size_t(n2))}, nil

}

func New(n1 int, n2 int) (*Matrix, error) {
	if n1 < 1 || n2 < 1  {
		return nil, errors.New("matrix dimensions should be > 0")
	}
	data := make([]float64, n1 * n2)
	return &Matrix{data, C.gsl_matrix_view_array((*C.double)(&data[0]),
		C.size_t(n1), C.size_t(n2))}, nil
}

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
			buffer.WriteString(strconv.FormatFloat(m.data[int(i1) * int(m.matrixView.matrix.size2) + int(i2)], 'E', 3, 64))
			if i2 < m.matrixView.matrix.size2 - 1 {
				buffer.WriteByte('\t')
			}
		}
		buffer.WriteByte('\n')
	}
	buffer.WriteByte('>')
	return buffer.String()
}

func (m *Matrix) Empty() *Matrix {
	nm, _ := New(		
		int(m.matrixView.matrix.size1),
		int(m.matrixView.matrix.size2))
	return nm
}

func (m *Matrix) Copy(dest *Matrix) {
	C.gsl_matrix_memcpy(&dest.matrixView.matrix, &m.matrixView.matrix)
}
func (m *Matrix) Scale(x float64) {
	C.gsl_matrix_scale(&m.matrixView.matrix, C.double(x))
}

func (m *Matrix) Exponential(em *Matrix) {
	C.gsl_linalg_exponential_ss(&m.matrixView.matrix, &em.matrixView.matrix, 0);
}

func (m *Matrix) SetItem(i1, i2 int, x float64) {
	m.data[i1 * int(m.matrixView.matrix.size2) + i2] = x
}

func (m *Matrix) ScaleItem(i1, i2 int, x float64) {
	m.data[i1 * int(m.matrixView.matrix.size2) + i2] = m.data[i1 * int(m.matrixView.matrix.size2) + i2] * x
}

func (m *Matrix) GetItem(i1, i2 int) float64 {
	return m.data[i1 * int(m.matrixView.matrix.size2) + i2]
}

func (m *Matrix) GetSize() (int, int) {
	return int(m.matrixView.matrix.size1), int(m.matrixView.matrix.size2)
}

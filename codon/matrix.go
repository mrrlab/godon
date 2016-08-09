package codon

import (
	"fmt"
	"math"
	"sort"

	"github.com/gonum/matrix/mat64"

	"bitbucket.org/Davydov/godon/bio"
)

// smallScale is a small value such that if Q-scale is less than it,
// the matrix is replace by an identity matrix.
const smallScale = 1e-30

var (
	// ZeroQ is a zero-matrix.
	ZeroQ     *mat64.Dense
	// IdentityP is an identity matrix. It's used if the branch or
	// scale is very small.
	IdentityP *mat64.Dense
)

// createIdentityMatrix creates an identity matrix of size size.
func createIdentityMatrix(size int) (m *mat64.Dense) {
	m = mat64.NewDense(size, size, nil)
	for i := 0; i < size; i++ {
		m.Set(i, i, 1)
	}
	return
}

// Sum calculates matrix sum.
func Sum(m *mat64.Dense) (s float64) {
	rows, cols := m.Dims()
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			s += math.Abs(m.At(i, j))
		}
	}
	return
}

// CreateTransitionMatrix creates a transition matrix.
func CreateTransitionMatrix(cf CodonFrequency, kappa, omega float64, m *mat64.Dense) (*mat64.Dense, float64) {
	return CreateRateTransitionMatrix(cf, kappa, omega, []float64{1, 1, 1}, m)
}

// CreateRateTransitionMatrix creates a transition matrix given the vector of rates.
func CreateRateTransitionMatrix(cf CodonFrequency, kappa, omega float64, rates []float64, m *mat64.Dense) (*mat64.Dense, float64) {
	//fmt.Println("kappa=", kappa, ", omega=", omega)
	if m == nil {
		m = mat64.NewDense(NCodon, NCodon, nil)
	}
	for i1 := 0; i1 < NCodon; i1++ {
		for i2 := 0; i2 < NCodon; i2++ {
			if i1 == i2 {
				m.Set(i1, i2, 0)
				continue
			}
			c1 := NumCodon[byte(i1)]
			c2 := NumCodon[byte(i2)]
			dist, transitions, pos := codonDistance(c1, c2)

			if dist > 1 {
				m.Set(i1, i2, 0)
				continue
			}
			m.Set(i1, i2, rates[pos])
			m.Set(i1, i2, m.At(i1, i2)*cf[i2])
			if transitions == 1 {
				m.Set(i1, i2, m.At(i1, i2)*kappa)
			}
			if bio.GeneticCode[c1] != bio.GeneticCode[c2] {
				m.Set(i1, i2, m.At(i1, i2)*omega)
			}
		}
	}
	for i1 := 0; i1 < NCodon; i1++ {
		rowSum := 0.0
		for i2 := 0; i2 < NCodon; i2++ {
			rowSum += m.At(i1, i2)
		}
		m.Set(i1, i1, -rowSum)
	}
	scale := 0.0
	for i := 0; i < NCodon; i++ {
		scale += -cf[i] * m.At(i, i)
	}

	if scale < smallScale {
		return ZeroQ, 0
	}
	return m, scale

}

// scaleMatrix scles a matrix by multiplying every value by scale.
func scaleMatrix(in *mat64.Dense, scale float64, out *mat64.Dense) *mat64.Dense {
	if out == nil {
		out = mat64.NewDense(NCodon, NCodon, nil)
	}
	out.Scale(scale, in)
	return out
}

// PrintQ prints a Q or P matrix.
func PrintQ(Q *mat64.Dense) {
	codons := make([]string, len(CodonNum))
	i := 0
	for k, _ := range CodonNum {
		codons[i] = k
		i++
	}
	sort.Strings(codons)

	fmt.Print("\t")
	for _, codon := range codons {
		fmt.Print(codon, "\t")
	}
	fmt.Println()
	for _, codon1 := range codons {
		fmt.Print(codon1, "\t")
		for _, codon2 := range codons {
			fmt.Printf("%0.4f\t", Q.At(int(CodonNum[codon1]), int(CodonNum[codon2])))
		}
		fmt.Println()
	}
}

// PrintUnQ prints Q or P matrix without codon names.
func PrintUnQ(Q *mat64.Dense) {
	fmt.Print("\t")
	for i := 0; i < NCodon; i++ {
		fmt.Print(NumCodon[byte(i)], "\t")
	}
	fmt.Println()
	for i1 := 0; i1 < NCodon; i1++ {
		fmt.Print(NumCodon[byte(i1)], "\t")
		for i2 := 0; i2 < NCodon; i2++ {
			fmt.Printf("%0.4f\t", Q.At(i1, i2))
		}
		fmt.Println()
	}
}

// codonDistance computes distance, number of transitions and
// difference position for two codons.
func codonDistance(c1, c2 string) (dist, transitions, pos int) {
	pos = -1
	for i := 0; i < len(c1); i++ {
		s1 := c1[i]
		s2 := c2[i]
		if s1 != s2 {
			pos = i
			dist++
			if ((s1 == 'A' || s1 == 'G') && (s2 == 'A' || s2 == 'G')) ||
				((s1 == 'T' || s1 == 'C') && (s2 == 'T' || s2 == 'C')) {
				transitions++
			}
		}
	}
	return
}

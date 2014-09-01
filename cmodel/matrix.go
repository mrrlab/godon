package cmodel

import (
	"fmt"
	"math"
	"sort"

	"github.com/skelterjohn/go.matrix"

	"bitbucket.com/Davydov/golh/bio"
)

// Sum calculates matrix sum.
func Sum(m *matrix.DenseMatrix) (s float64) {
	for i := 0; i < m.Rows(); i++ {
		for j := 0; j < m.Cols(); j++ {
			s += math.Abs(m.Get(i, j))
		}
	}
	return
}

func createTransitionMatrix(cf CodonFrequency, kappa, omega float64, m *matrix.DenseMatrix) (*matrix.DenseMatrix, float64) {
	//fmt.Println("kappa=", kappa, ", omega=", omega)
	if m == nil {
		m = matrix.Zeros(nCodon, nCodon)
	}
	for i1 := 0; i1 < nCodon; i1++ {
		for i2 := 0; i2 < nCodon; i2++ {
			if i1 == i2 {
				m.Set(i1, i2, 0)
				continue
			}
			c1 := numCodon[byte(i1)]
			c2 := numCodon[byte(i2)]
			dist, transitions := codonDistance(c1, c2)

			if dist > 1 {
				m.Set(i1, i2, 0)
				continue
			}
			m.Set(i1, i2, cf[i2])
			if transitions == 1 {
				m.Set(i1, i2, m.Get(i1, i2)*kappa)
			}
			if bio.GeneticCode[c1] != bio.GeneticCode[c2] {
				m.Set(i1, i2, m.Get(i1, i2)*omega)
			}
		}
	}
	for i1 := 0; i1 < nCodon; i1++ {
		rowSum := 0.0
		for i2 := 0; i2 < nCodon; i2++ {
			rowSum += m.Get(i1, i2)
		}
		m.Set(i1, i1, -rowSum)
	}
	scale := 0.0
	for i := 0; i < nCodon; i++ {
		scale += -cf[i] * m.Get(i, i)
	}

	return m, scale

}

func PrintQ(Q *matrix.DenseMatrix) {
	codons := make([]string, len(codonNum))
	i := 0
	for k, _ := range codonNum {
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
			fmt.Printf("%0.4f\t", Q.Get(int(codonNum[codon1]), int(codonNum[codon2])))
		}
		fmt.Println()
	}
}

func PrintUnQ(Q *matrix.DenseMatrix) {
	fmt.Print("\t")
	for i := 0; i < nCodon; i++ {
		fmt.Print(numCodon[byte(i)], "\t")
	}
	fmt.Println()
	for i1 := 0; i1 < nCodon; i1++ {
		fmt.Print(numCodon[byte(i1)], "\t")
		for i2 := 0; i2 < nCodon; i2++ {
			fmt.Printf("%0.4f\t", Q.Get(i1, i2))
		}
		fmt.Println()
	}
}

func codonDistance(c1, c2 string) (dist, transitions int) {
	for i := 0; i < len(c1); i++ {
		s1 := c1[i]
		s2 := c2[i]
		if s1 != s2 {
			dist++
			if ((s1 == 'A' || s1 == 'G') && (s2 == 'A' || s2 == 'G')) ||
				((s1 == 'T' || s1 == 'C') && (s2 == 'T' || s2 == 'C')) {
				transitions++
			}
		}
	}
	return
}

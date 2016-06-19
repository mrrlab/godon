package codon

import (
	"github.com/gonum/matrix/mat64"
	"testing"
)

func BenchmarkEigen1(b *testing.B) {
	cf := F0()
	q := mat64.NewDense(NCodon, NCodon, nil)
	p := mat64.NewDense(NCodon, NCodon, nil)
	for i := 0; i < b.N; i++ {
		q, s := CreateTransitionMatrix(cf, 2.1, 0.25, q)
		e := EMatrix{Q: q, Scale: s, CF: cf}
		e.Eigen()
		e.Exp(p, 0.3)
	}
}

func BenchmarkEigen2(b *testing.B) {
	cf := F0()
	q := mat64.NewDense(NCodon, NCodon, nil)
	p := mat64.NewDense(NCodon, NCodon, nil)
	q, s := CreateTransitionMatrix(cf, 2.1, 0.25, q)
	e := EMatrix{Q: q, Scale: s, CF: cf}
	e.Eigen()
	for i := 0; i < b.N; i++ {
		e.Exp(p, 0.3)
	}
}

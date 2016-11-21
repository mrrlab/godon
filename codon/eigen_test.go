package codon

import (
	"github.com/gonum/matrix/mat64"
	"testing"

	"bitbucket.org/Davydov/godon/bio"
)

func BenchmarkEigen1(b *testing.B) {
	// this is a hack to pass genetic code to F0
	gcode := bio.GeneticCodes[1]
	NCodon := gcode.NCodon
	cs := []Sequence{
		Sequence{GCode: gcode},
	}

	cf := F0(cs)
	q := mat64.NewDense(NCodon, NCodon, nil)
	p := mat64.NewDense(NCodon, NCodon, nil)
	for i := 0; i < b.N; i++ {
		q, s := CreateTransitionMatrix(cf, 2.1, 0.25, q)
		e := NewEMatrix(cf)
		e.Set(q, s)
		e.Eigen()
		e.Exp(p, 0.3)
	}
}

func BenchmarkEigen2(b *testing.B) {
	gcode := bio.GeneticCodes[1]
	NCodon := gcode.NCodon
	cs := []Sequence{
		Sequence{GCode: gcode},
	}

	cf := F0(cs)
	q := mat64.NewDense(NCodon, NCodon, nil)
	p := mat64.NewDense(NCodon, NCodon, nil)
	q, s := CreateTransitionMatrix(cf, 2.1, 0.25, q)
	e := NewEMatrix(cf)
	e.Set(q, s)
	e.Eigen()
	for i := 0; i < b.N; i++ {
		e.Exp(p, 0.3)
	}
}

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
		{GCode: gcode},
	}

	cf := F0(cs)
	q := mat64.NewDense(NCodon, NCodon, nil)
	p := mat64.NewDense(NCodon, NCodon, nil)
	for i := 0; i < b.N; i++ {
		var s float64
		q, s = CreateTransitionMatrix(cf, 2.1, 0.25, q)
		e := NewEMatrix(cf)
		e.Set(q, s)
		if err := e.Eigen(); err != nil {
			b.Error("Error: ", err)
		}
		tmp := make([]float64, NCodon*NCodon)
		res := make([]float64, NCodon*NCodon)
		if _, err := e.Exp(p, 0.3, res, tmp); err != nil {
			b.Error("Error: ", err)
		}
	}
}

func BenchmarkEigen2(b *testing.B) {
	gcode := bio.GeneticCodes[1]
	NCodon := gcode.NCodon
	cs := []Sequence{
		{GCode: gcode},
	}

	cf := F0(cs)
	q := mat64.NewDense(NCodon, NCodon, nil)
	p := mat64.NewDense(NCodon, NCodon, nil)
	q, s := CreateTransitionMatrix(cf, 2.1, 0.25, q)
	e := NewEMatrix(cf)
	e.Set(q, s)
	if err := e.Eigen(); err != nil {
		b.Error("Error: ", err)
	}
	for i := 0; i < b.N; i++ {
		tmp := make([]float64, NCodon*NCodon)
		res := make([]float64, NCodon*NCodon)
		if _, err := e.Exp(p, 0.3, res, tmp); err != nil {
			b.Error("Error: ", err)
		}
	}
}

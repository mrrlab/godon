package paml

import (
	"math"
	"testing"
)

const smallDiff = 1e-6

type Settings struct {
	n      int
	a, b   float64
	median bool
}

/*** Tests if a and b are approximately equal ***/
func appreq(a, b float64) bool {
	if math.Abs(a-b) > smallDiff {
		return false
	}
	return true
}

/*** Tests that arrays have approximately same values ***/
func cmp(a, b []float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !appreq(a[i], b[i]) {
			return false
		}
	}
	return true
}

/*** Tests that all values are equal and sum to 1 ***/
func alleq(r []float64) bool {
	if len(r) < 1 {
		return false
	}
	v := r[0]
	sum := v
	for i := 1; i < len(r); i++ {
		if !appreq(v, r[i]) {
			return false
		}
		sum += r[i]
	}
	if !appreq(sum, 1) {
		return false
	}
	return true
}

/*** Test discrete beta ***/
func TestBeta(tst *testing.T) {
	settings := [...]Settings{
		Settings{4, 0.5, 10, false},
		Settings{4, 0.5, 10, true},
		Settings{8, 2, .1, false},
		Settings{7, 15, 1, true},
	}
	results := [...]([]float64){
		[]float64{0.001709, 0.012818, 0.041099, 0.134851},
		[]float64{0.001449, 0.013916, 0.045204, 0.129907},
		[]float64{0.687057, 0.943958, 0.989629, 0.998540, 0.999869, 0.999994, 1.000000, 1.000000},
		[]float64{0.836482, 0.900046, 0.931225, 0.952350, 0.968440, 0.981483, 0.992475},
	}
	for i, s := range settings {
		freq := make([]float64, s.n)
		r := DiscreteBeta(s.a, s.b, s.n, s.median, freq, nil)
		if !alleq(freq) {
			tst.Error("Frequencies are wrong:", freq)
		}
		if !cmp(r, results[i]) {
			tst.Error("Results missmatch:", r, results[i])
		}
	}
}

/*** Test discrete gamma ***/
func TestGamma(tst *testing.T) {
	settings := [...]Settings{
		Settings{4, 0.5, 10, false},
		Settings{4, 0.5, 10, true},
		Settings{8, 2, .1, false},
		Settings{7, 15, 1, true},
	}
	results := [...]([]float64){
		[]float64{0.001669, 0.012596, 0.041013, 0.144721},
		[]float64{0.001454, 0.014036, 0.046239, 0.138272},
		[]float64{3.848344, 7.882645, 11.320993, 14.879554, 18.906079, 23.893507, 31.028044, 48.240834},
		[]float64{9.793787, 11.891047, 13.362596, 14.722906, 16.172736, 17.973174, 21.083754},
	}
	for i, s := range settings {
		freq := make([]float64, s.n)
		r := DiscreteGamma(s.a, s.b, s.n, s.median, freq, nil)
		if !alleq(freq) {
			tst.Error("Frequencies are wrong:", freq)
		}
		if !cmp(r, results[i]) {
			tst.Error("Results missmatch:", r, results[i])
		}
	}
}

package optimize

import (
	"math/rand"
)

// Returns a random value in the range [0, 1], including 1.
func Rand() float64 {
	// 1.0 is not included and we would like to be symmetric
	r := float64(1)
	for r > 0.999 {
		r = rand.Float64()
	}
	return r / 0.999

}

// UniformProposal returns uniform proposal function.
func UniformProposal(width float64) func(float64) float64 {
	if width <= 0 {
		panic("width should be non-negative")
	}
	return func(x float64) float64 {
		return x + Rand()*width - width/2
	}
}

// UniformGlobalProposal returns uniform proposal function given max
// and min.
func UniformGlobalProposal(min, max float64) func(float64) float64 {
	if max <= min {
		panic("max <= min")
	}
	return func(x float64) float64 {
		return Rand()*(max-min) - min
	}
}

// NormalProposal returns normal proposal function.
func NormalProposal(sd float64) func(float64) float64 {
	if sd <= 0 {
		panic("sd should be >= 0")
	}
	return func(x float64) float64 {
		return x + rand.NormFloat64()*sd
	}
}

// DiscreteProposal returns function returning a random integer
// converted to float64.
func DiscreteProposal(state int, nstates int) (newstate int) {
	if nstates <= 1 {
		panic("number of states should be at least 1")
	}
	if state < 0 {
		panic("incorrect state")
	}
	newstate = rand.Intn(nstates - 1)
	if newstate >= state {
		newstate++
	}
	return
}

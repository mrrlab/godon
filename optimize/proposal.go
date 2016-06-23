package optimize

import (
	"math/rand"
)

func Rand() float64 {
	// 1.0 is not included and we would like to be symmetric
	r := float64(1)
	for r > 0.999 {
		r = rand.Float64()
	}
	return r / 0.999

}

func UniformProposal(width float64) func(float64) float64 {
	if width <= 0 {
		panic("width should be non-negative")
	}
	return func(x float64) float64 {
		return x + Rand()*width - width/2
	}
}

func UniformGlobalProposal(min, max float64) func(float64) float64 {
	if max <= min {
		panic("max <= min")
	}
	return func(x float64) float64 {
		return Rand()*(max-min) - min
	}
}

func NormalProposal(sd float64) func(float64) float64 {
	if sd <= 0 {
		panic("sd should be >= 0")
	}
	return func(x float64) float64 {
		return x + rand.NormFloat64()*sd
	}
}

func DiscretePropose(state int, nstates int) (newstate int) {
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

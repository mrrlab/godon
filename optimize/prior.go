package optimize

import (
	"math"
)

func UniformPrior(min, max float64, incmin, incmax bool) func(float64) float64 {
	if max <= min {
		panic("max <= min")
	}
	return func(x float64) float64 {
		if (incmin && x < min) ||
			(!incmin && x <= min) ||
			(incmax && x > max) ||
			(!incmax && x >= max) {
			return math.Inf(-1)
		} else {
			return -math.Log(max - min)
		}
	}
}

func GammaPrior(shape, scale float64, inczero bool) func(float64) float64 {
	if scale <= 0 || scale <= 0 {
		panic("shape and scale of gamma distribution must be > 0")
	}
	return func(x float64) float64 {
		if x < 0 || (x == 0 && !inczero) {
			return math.Inf(-1)
		}
		g, _ := math.Lgamma(shape)
		return (shape-1)*math.Log(x) - x/scale - shape*math.Log(scale) - g
	}
}

func ExponentialPrior(rate float64, inczero bool) func(float64) float64 {
	if rate <= 0 {
		panic("exponential rate should be > 0")
	}
	return func(x float64) float64 {
		if x < 0 || (x == 0 && !inczero) {
			return math.Inf(-1)
		}
		return math.Log(rate) - rate*x
	}
}

func ProductPrior(f, g func(float64) float64) func(float64) float64 {
	return func(x float64) float64 {
		return f(x) * g(x)
	}
}

package mcmc

import (
	"fmt"
	"math"
	"math/rand"
)

type Parameter interface {
	Prior() float64
	Propose()
	Reject()
	Value() fmt.Stringer
}

type Float64Parameter struct {
	float64
	old          float64
	PriorFunc    func(float64) float64
	ProposalFunc func(float64) float64
	Min          float64
	Max          float64
}

func (p *Float64Parameter) Prior() float64 {
	return p.PriorFunc(p.float64)
}

func (p *Float64Parameter) Propose() {
	p.old, p.float64 = p.float64, p.ProposalFunc(p.float64)

	for p.float64 < p.Min || p.float64 > p.Max {
		if p.float64 < p.Min {
			p.float64 = p.Min + (p.Min - p.float64)
		}
		if p.float64 > p.Max {
			p.float64 = p.Max - (p.float64 - p.Max)
		}
	}
}

func (p *Float64Parameter) Reject() {
	p.float64 = p.old
}

func (p *Float64Parameter) Value() float64 {
	return p.float64
}

func FlatPrior(min, max float64, incmin, incmax bool) func(float64) float64 {
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

func ProductPrior(f, g func(float64) float64) func(float64) float64 {
	return func(x float64) float64 {
		return f(x) * g(x)
	}
}

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

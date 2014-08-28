package mcmc

import (
	"math"
	"math/rand"
	"strconv"
)

type Parameter interface {
	Name() string
	Prior() float64
	OldPrior() float64
	Propose()
	Accept()
	Reject()
	Value() string
}

type Parameters []Parameter

type Float64Parameters []*Float64Parameter

func (p Float64Parameters) Get(n string) float64 {
	for _, par := range p {
		if par.Name() == n {
			return *par.float64
		}
	}
	panic("unknown parameter name")
	//return math.NaN()
}

func (p Float64Parameters) Set(n string, x float64) {
	for _, par := range p {
		if par.Name() == n {
			par.old, *par.float64 = *par.float64, x
			par.reflect()
			par.OnChange()
			return
		}
	}
	panic("unknown parameter name")
}

type Float64Parameter struct {
	*float64
	old          float64
	name         string
	PriorFunc    func(float64) float64
	ProposalFunc func(float64) float64
	Min          float64
	Max          float64
	OnChange     func()
}

func NewFloat64Parameter(par *float64, name string) *Float64Parameter {
	return &Float64Parameter{
		float64:      par,
		name:         name,
		PriorFunc:    UniformPrior(-1, 1, true, true),
		ProposalFunc: NormalProposal(1),
		Min:          math.Inf(-1),
		Max:          math.Inf(+1),
	}
}

func (p *Float64Parameter) Name() string {
	return p.name
}

func (p *Float64Parameter) Prior() float64 {
	return p.PriorFunc(*p.float64)
}

func (p *Float64Parameter) OldPrior() float64 {
	return p.PriorFunc(p.old)
}

func (p *Float64Parameter) reflect() {
	for *p.float64 < p.Min || *p.float64 > p.Max {
		if *p.float64 < p.Min {
			*p.float64 = p.Min + (p.Min - *p.float64)
		}
		if *p.float64 > p.Max {
			*p.float64 = p.Max - (*p.float64 - p.Max)
		}
	}
}

func (p *Float64Parameter) Propose() {
	p.old, *p.float64 = *p.float64, p.ProposalFunc(*p.float64)
	p.reflect()
	if p.OnChange != nil {
		p.OnChange()
	}
}

func (p *Float64Parameter) Reject() {
	*p.float64, p.old = p.old, *p.float64
	if p.OnChange != nil {
		p.OnChange()
	}
}

func (p *Float64Parameter) Accept() {
}

func (p *Float64Parameter) Value() string {
	return strconv.FormatFloat(*p.float64, 'f', 6, 64)
}

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
		if x < 0 || x == 0 || !inczero {
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

package optimize

import (
	"fmt"
	"math"
	"strconv"
)

type Parameter interface {
	Name() string
	Prior() float64
	OldPrior() float64
	Propose()
	Accept(int)
	Reject()
	GetValue() fmt.Stringer
	SetValue(interface{})
}

type Parameters []Parameter

type Float64Parameters []*Float64Parameter

type Float64 float64

func (f *Float64) String() string {
	return strconv.FormatFloat(*(*float64)(f), 'f', 6, 64)
}

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

func (p *Float64Parameter) Accept(iter int) {
}

func (p *Float64Parameter) GetValue() fmt.Stringer {
	return (*Float64)(p.float64)
}

func (p *Float64Parameter) SetValue(x interface{}) {
	*p.float64 = x.(float64)
}

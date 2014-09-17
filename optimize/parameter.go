package optimize

import (
	"math"
	"strconv"
)

type FloatParameter interface {
	Name() string
	Prior() float64
	OldPrior() float64
	Propose()
	Accept(int)
	Reject()
	String() string
	Get() float64
	Set(float64)
}

type FloatParameters []FloatParameter

func (p *FloatParameters) Append(par FloatParameter) {
	*p = append(*p, par)
}

func (p *FloatParameters) Names(is []string) (s []string) {
	if is == nil {
		s = make([]string, len(*p))
	} else {
		s = is
	}
	for i, par := range *p {
		s[i] = par.Name()
	}
	return
}

func (p *FloatParameters) Values(iv []float64) (v []float64) {
	if v == nil {
		v = make([]float64, len(*p))
	} else {
		v = iv
	}
	for i, par := range *p {
		v[i] = par.Get()
	}
	return
}

func (p *FloatParameters) NamesString() (s string) {
	for i, par := range *p {
		if i != 0 {
			s += "\t"
		}
		s += par.Name()
	}
	return
}

func (p *FloatParameters) ValuesString() (s string) {
	for i, par := range *p {
		if i != 0 {
			s += "\t"
		}
		s += par.String()
	}
	return
}

type BasicFloatParameter struct {
	*float64
	old          float64
	name         string
	PriorFunc    func(float64) float64
	ProposalFunc func(float64) float64
	Min          float64
	Max          float64
	OnChange     func()
}

func NewBasicFloatParameter(par *float64, name string) *BasicFloatParameter {
	return &BasicFloatParameter{
		float64:      par,
		name:         name,
		PriorFunc:    UniformPrior(-1, 1, true, true),
		ProposalFunc: NormalProposal(1),
		Min:          math.Inf(-1),
		Max:          math.Inf(+1),
	}
}

func (p *BasicFloatParameter) Get() float64 {
	return *p.float64
}

func (p *BasicFloatParameter) Set(v float64) {
	*p.float64 = v
	if p.OnChange != nil {
		p.OnChange()
	}
}
func (p *BasicFloatParameter) Name() string {
	return p.name
}

func (p *BasicFloatParameter) Prior() float64 {
	return p.PriorFunc(*p.float64)
}

func (p *BasicFloatParameter) OldPrior() float64 {
	return p.PriorFunc(p.old)
}

func (p *BasicFloatParameter) reflect() {
	for *p.float64 < p.Min || *p.float64 > p.Max {
		if *p.float64 < p.Min {
			*p.float64 = p.Min + (p.Min - *p.float64)
		}
		if *p.float64 > p.Max {
			*p.float64 = p.Max - (*p.float64 - p.Max)
		}
	}
}

func (p *BasicFloatParameter) Propose() {
	p.old, *p.float64 = *p.float64, p.ProposalFunc(*p.float64)
	p.reflect()
	if p.OnChange != nil {
		p.OnChange()
	}
}

func (p *BasicFloatParameter) Reject() {
	*p.float64, p.old = p.old, *p.float64
	if p.OnChange != nil {
		p.OnChange()
	}
}

func (p *BasicFloatParameter) Accept(iter int) {
}

func (p *BasicFloatParameter) String() string {
	return strconv.FormatFloat(*p.float64, 'f', 6, 64)
}

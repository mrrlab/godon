package optimize

import (
	"errors"
	"math"
	"math/rand"
	"strconv"
)

const (
	MIN = -10
	MAX = +10
)

type FloatParameter interface {
	Name() string
	Prior() float64
	OldPrior() float64
	Propose()
	Accept(int)
	Reject()
	String() string
	SetMin(float64)
	SetMax(float64)
	GetMin() float64
	GetMax() float64
	SetOnChange(func())
	SetProposalFunc(func(float64) float64)
	SetPriorFunc(func(float64) float64)
	Get() float64
	Set(float64)
	InRange() bool
	ValueInRange(float64) bool
}

type FloatParameterGenerator func(*float64, string) FloatParameter

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

func (p *FloatParameters) ValuesInRange(vals []float64) bool {
	if len(vals) != len(*p) {
		panic("Incorrect number of parameters")
	}
	for i, par := range *p {
		if !par.ValueInRange(vals[i]) {
			return false
		}
	}
	return true
}

func (p *FloatParameters) SetValues(v []float64) error {
	if len(v) != len(*p) {
		errors.New("Incorrect number of parameters")
	}
	for i, par := range *p {
		par.Set(v[i])
	}
	return nil
}

func (par *FloatParameters) ReadLine(l string) error {
	v, err := ReadFloats(l)
	if err != nil {
		return err
	}
	err = par.SetValues(v[2:])
	if err != nil {
		return err
	}
	return nil
}

func (p *FloatParameters) Update(pSrc *FloatParameters) {
	for i := range *p {
		(*p)[i].Set((*pSrc)[i].Get())
	}
}

func (p *FloatParameters) Randomize() {
	for _, par := range *p {
		min := math.Max(MIN, par.GetMin())
		max := math.Min(MAX, par.GetMax())
		d := max - min
		par.Set(min + rand.Float64()*d)
	}

}
func (p *FloatParameters) InRange() bool {
	for _, par := range *p {
		if !par.InRange() {
			return false
		}
	}
	return true
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
	priorFunc    func(float64) float64
	proposalFunc func(float64) float64
	min          float64
	max          float64
	onChange     func()
}

func NewBasicFloatParameter(par *float64, name string) *BasicFloatParameter {
	return &BasicFloatParameter{
		float64:      par,
		name:         name,
		priorFunc:    UniformPrior(-1, 1, true, true),
		proposalFunc: NormalProposal(1),
		min:          math.Inf(-1),
		max:          math.Inf(+1),
	}
}

func BasicFloatParameterGenerator(par *float64, name string) FloatParameter {
	return NewBasicFloatParameter(par, name)
}

func (p *BasicFloatParameter) SetMin(min float64) {
	p.min = min
}

func (p *BasicFloatParameter) SetMax(max float64) {
	p.max = max
}

func (p *BasicFloatParameter) SetPriorFunc(f func(float64) float64) {
	p.priorFunc = f
}

func (p *BasicFloatParameter) SetProposalFunc(f func(float64) float64) {
	p.proposalFunc = f
}

func (p *BasicFloatParameter) SetOnChange(f func()) {
	p.onChange = f
}

func (p *BasicFloatParameter) Get() float64 {
	return *p.float64
}

func (p *BasicFloatParameter) Set(v float64) {
	if *p.float64 == v {
		// do nothing if value has not cnahged
		return
	}
	*p.float64 = v
	if p.onChange != nil {
		p.onChange()
	}
}

func (p *BasicFloatParameter) GetMin() float64 {
	return p.min
}

func (p *BasicFloatParameter) GetMax() float64 {
	return p.max
}

func (p *BasicFloatParameter) ValueInRange(v float64) bool {
	if v < p.min || v > p.max {
		return false
	}
	return true
}

func (p *BasicFloatParameter) InRange() bool {
	if *p.float64 < p.min || *p.float64 > p.max {
		return false
	}
	return true
}

func (p *BasicFloatParameter) Name() string {
	return p.name
}

func (p *BasicFloatParameter) Prior() float64 {
	return p.priorFunc(*p.float64)
}

func (p *BasicFloatParameter) OldPrior() float64 {
	return p.priorFunc(p.old)
}

func (p *BasicFloatParameter) reflect() {
	for *p.float64 < p.min || *p.float64 > p.max {
		if *p.float64 < p.min {
			*p.float64 = p.min + (p.min - *p.float64)
		}
		if *p.float64 > p.max {
			*p.float64 = p.max - (*p.float64 - p.max)
		}
	}
}

func (p *BasicFloatParameter) Propose() {
	p.old, *p.float64 = *p.float64, p.proposalFunc(*p.float64)
	p.reflect()
	if p.onChange != nil {
		p.onChange()
	}
}

func (p *BasicFloatParameter) Reject() {
	*p.float64, p.old = p.old, *p.float64
	if p.onChange != nil {
		p.onChange()
	}
}

func (p *BasicFloatParameter) Accept(iter int) {
}

func (p *BasicFloatParameter) String() string {
	return strconv.FormatFloat(*p.float64, 'f', 6, 64)
}

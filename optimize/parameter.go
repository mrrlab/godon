package optimize

import (
	"errors"
	"math"
	"math/rand"
	"strconv"
)

const (
	// Minimum value for randomization.
	MIN = -10
	// Maximum value for randomization.
	MAX = +10
)

// FloatParameter is a floating point parameter.
type FloatParameter interface {
	// Name returns the parameter name.
	Name() string
	// Prior returns prior value for the current x.
	Prior() float64
	// OldPrior returns prior value for the old x.
	OldPrior() float64
	// Propose set x to a new proposed value
	Propose()
	// Accept accepts the proposal.
	Accept(int)
	// Reject rejects the proposal (and set x to an old value).
	Reject()
	// String returns a string representation of a value.
	String() string
	// SetMin sets a minimal value.
	SetMin(float64)
	// SetMax sets a maximal value.
	SetMax(float64)
	// GetMin returns a minimal value.
	GetMin() float64
	// GetMax returns a maximal value.
	GetMax() float64
	// SetOnChange sets a callback which should be called when the
	// value is changed.
	SetOnChange(func())
	// SetProposalFunc sets a proposal for the parameter.
	SetProposalFunc(func(float64) float64)
	// SetPriorFunc sets a prior for the parameter.
	SetPriorFunc(func(float64) float64)
	// Get returns the parameter value.
	Get() float64
	// Set sets the parameter.
	Set(float64)
	// InRange returns true if parameter value is between min and
	// max.
	InRange() bool
	// ValueInRange returns true if value is between min and max.
	ValueInRange(float64) bool
}

// FloatParameterGenerator is a function which returns a
// FloatParameter.
type FloatParameterGenerator func(*float64, string) FloatParameter

// FloatParameters is an array of FloatParameters.
type FloatParameters []FloatParameter

// Append appends a new float parameter.
func (p *FloatParameters) Append(par FloatParameter) {
	*p = append(*p, par)
}

// Names returns a slice of parameter names.
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

// Values returns a slice of parameter values.
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

// ValuesInRange returns true if all parameters are in range.
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

// SetValues sets values from a slice.
func (p *FloatParameters) SetValues(v []float64) error {
	if len(v) != len(*p) {
		errors.New("Incorrect number of parameters")
	}
	for i, par := range *p {
		par.Set(v[i])
	}
	return nil
}

// ReadLine sets values from a string.
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

// Update updates values from another FloatParameters object.
func (p *FloatParameters) Update(pSrc *FloatParameters) {
	for i := range *p {
		(*p)[i].Set((*pSrc)[i].Get())
	}
}

// Randomize sets random values to all the parameters. If no minimum
// and maximum are specified, MIN and MAX constants are used.
func (p *FloatParameters) Randomize() {
	for _, par := range *p {
		min := math.Max(MIN, par.GetMin())
		max := math.Min(MAX, par.GetMax())
		d := max - min
		par.Set(min + rand.Float64()*d)
	}

}

// InRange returns true if parameter value is between min and
// max.
func (p *FloatParameters) InRange() bool {
	for _, par := range *p {
		if !par.InRange() {
			return false
		}
	}
	return true
}

// Returns a string with tab-separated parameter names.
func (p *FloatParameters) NamesString() (s string) {
	for i, par := range *p {
		if i != 0 {
			s += "\t"
		}
		s += par.Name()
	}
	return
}

// ValuesString returns a tab-separated string of values.
func (p *FloatParameters) ValuesString() (s string) {
	for i, par := range *p {
		if i != 0 {
			s += "\t"
		}
		s += par.String()
	}
	return
}

// BasicFloatParameter implementation of FloatParameter interface.
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

// NewBasicFloatParameter creates a new BasicFloatParameter.
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

// BasicFloatParameterGenerator is a function which returns a
// FloatParameter which is a BasicFloatParameter.
func BasicFloatParameterGenerator(par *float64, name string) FloatParameter {
	return NewBasicFloatParameter(par, name)
}

// SetMin sets a minimal value.
func (p *BasicFloatParameter) SetMin(min float64) {
	p.min = min
}

// SetMax sets a maximal value.
func (p *BasicFloatParameter) SetMax(max float64) {
	p.max = max
}

// SetPriorFunc sets a prior for the parameter.
func (p *BasicFloatParameter) SetPriorFunc(f func(float64) float64) {
	p.priorFunc = f
}

// SetProposalFunc sets a proposal for the parameter.
func (p *BasicFloatParameter) SetProposalFunc(f func(float64) float64) {
	p.proposalFunc = f
}

// SetOnChange sets a callback which should be called when the
// value is changed.
func (p *BasicFloatParameter) SetOnChange(f func()) {
	p.onChange = f
}

// Get returns the parameter value.
func (p *BasicFloatParameter) Get() float64 {
	return *p.float64
}

// Set sets the parameter.
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

// GetMin returns a minimal value.
func (p *BasicFloatParameter) GetMin() float64 {
	return p.min
}

// GetMax returns a maximum value.
func (p *BasicFloatParameter) GetMax() float64 {
	return p.max
}

// ValueInRange returns true if value is between min and max.
func (p *BasicFloatParameter) ValueInRange(v float64) bool {
	if v < p.min || v > p.max {
		return false
	}
	return true
}

// InRange returns true if parameter value is between min and
// max.
func (p *BasicFloatParameter) InRange() bool {
	if *p.float64 < p.min || *p.float64 > p.max {
		return false
	}
	return true
}

// Name returns the parameter name.
func (p *BasicFloatParameter) Name() string {
	return p.name
}

// Prior returns prior value for the current x.
func (p *BasicFloatParameter) Prior() float64 {
	return p.priorFunc(*p.float64)
}

// OldPrior returns prior value for the old x.
func (p *BasicFloatParameter) OldPrior() float64 {
	return p.priorFunc(p.old)
}

// reflect puts value into the range by iteratively reflecting from
// the borders.
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

// Propose set x to a new proposed value
func (p *BasicFloatParameter) Propose() {
	p.old, *p.float64 = *p.float64, p.proposalFunc(*p.float64)
	p.reflect()
	if p.onChange != nil {
		p.onChange()
	}
}

// Reject rejects the proposal (and set x to an old value).
func (p *BasicFloatParameter) Reject() {
	*p.float64, p.old = p.old, *p.float64
	if p.onChange != nil {
		p.onChange()
	}
}

// Accept accepts the proposal.
func (p *BasicFloatParameter) Accept(iter int) {
}

// String returns a string representation of a value.
func (p *BasicFloatParameter) String() string {
	return strconv.FormatFloat(*p.float64, 'f', 6, 64)
}

// DiscreteParameter is a discrete parameter. It is internally based
// on a BasicFloatParameter. The state is stored as an integer number
// stored converted to float (e.g. 3.0 instead of 3).
type DiscreteParameter struct {
	*BasicFloatParameter
	// NStates number of discrete states of the parameter.
	NStates int
}

// Creates a new DiscreteParameter.
func NewDiscreteParameter(par *float64, name string, nstates int) (dpar *DiscreteParameter) {
	if nstates <= 1 {
		panic("incorrect number of states for a discrete parameter")
	}
	dpar = &DiscreteParameter{NStates: nstates}
	dpar.BasicFloatParameter = &BasicFloatParameter{
		float64:   par,
		name:      name,
		priorFunc: DiscreteUniformPrior(nstates),
		min:       0,
		max:       float64(nstates - 1),
	}
	return
}

// Propose set x to a new proposed value
func (p *DiscreteParameter) Propose() {
	newState := DiscretePropose(int(*p.float64), p.NStates)
	p.old, *p.float64 = *p.float64, float64(newState)
	if p.onChange != nil {
		p.onChange()
	}
}

// String returns a string representation of a value.
func (p *DiscreteParameter) String() string {
	return strconv.Itoa(int(*p.float64))
}

package optimize

// None is an optimizer which computes initial value and exits.
type None struct {
	BaseOptimizer
}

// NewNone creates an optimizer which computes initial likelihood only.
func NewNone() *None {
	return &None{}
}

// Run starts sampling.
func (i *None) Run(iterations int) {
	i.l = i.Likelihood()
	i.calls++
	i.maxL = i.l
	i.maxLPar = i.parameters.Values(i.maxLPar)
	i.PrintHeader()
	i.PrintLine(i.parameters, i.l)
}

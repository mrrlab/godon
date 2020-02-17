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
	i.SaveStart()
	i.PrintHeader()
	i.PrintLine(i.parameters, i.startL, i.repPeriod)
	i.SaveCheckpoint(true)
	i.saveDeltaT()
}

package optimize

import (
	"math"

	lbfgsb "github.com/afbarnard/go-lbfgsb"
)

type LBFGSB struct {
	BaseOptimizer
	dH             float64
	grad           []float64
	suppressHeader bool
}

func NewLBFGSB() (lbfgsb *LBFGSB) {
	lbfgsb = &LBFGSB{
		BaseOptimizer: BaseOptimizer{
			repPeriod: 10,
		},
		dH: 1e-6,
	}
	return
}

func (l *LBFGSB) Logger(info *lbfgsb.OptimizationIterationInformation) {
	l.i += 1
	l.parameters.SetValues(info.X)
	l.PrintLine(l.parameters, -info.F)
	select {
	case s := <-l.sig:
		log.Fatal("Received signal exiting:", s)
	default:
	}
}

func (l *LBFGSB) EvaluateFunction(x []float64) float64 {
	if !l.parameters.ValuesInRange(x) {
		return math.Inf(+1)
	}

	l.parameters.SetValues(x)

	L := l.Likelihood()
	l.calls++
	if L > l.maxL {
		l.maxL = L
		l.maxLPar = l.parameters.Values(l.maxLPar)
	}
	return -L
}

func (l *LBFGSB) EvaluateGradient(x []float64) (grad []float64) {
	if l.grad == nil {
		l.grad = make([]float64, len(x))
	}
	grad = l.grad
	// we assume that values are in range
	l.parameters.SetValues(x)

	l1 := -l.Likelihood()
	l.calls++
	for i, _ := range x {
		inv := false
		v := x[i] + l.dH

		// this shouldn't happen with current boundaries
		// but to be safe
		if v >= l.parameters[i].GetMax() {
			v = x[i] - l.dH
			inv = true
		}

		l.parameters[i].Set(v)
		l2 := -l.Likelihood()
		l.calls++

		grad[i] = (l2 - l1) / l.dH
		if inv {
			grad[i] = -grad[i]
		}

		l.parameters[i].Set(x[i])

		select {
		case s := <-l.sig:
			log.Fatal("Received signal exiting:", s)
		default:
		}

	}
	return
}

func (l *LBFGSB) Run(iterations int) {
	l.maxL = math.Inf(-1)
	if !l.suppressHeader {
		l.PrintHeader(l.parameters)
	}
	bounds := make([][2]float64, len(l.parameters))

	for i, par := range l.parameters {
		bounds[i][0] = par.GetMin() + 1e-5
		bounds[i][1] = par.GetMax() - 1e-5
	}

	opt := new(lbfgsb.Lbfgsb)
	opt.SetApproximationSize(10)
	opt.SetFTolerance(1e-9)
	opt.SetGTolerance(1e-9)

	opt.SetBounds(bounds)
	opt.SetLogger(l.Logger)

	_, exitStatus := opt.Minimize(l, l.parameters.Values(nil))

	switch exitStatus.Code {
	case lbfgsb.SUCCESS:
		fallthrough
	case lbfgsb.APPROXIMATE:
		log.Info("LBFGSB status:", exitStatus)
	case lbfgsb.WARNING:
		log.Warning("Warning, LBFGSB status:", exitStatus)
	default:
		log.Error("Error during LBFGSB:", exitStatus)
	}
}

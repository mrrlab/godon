package optimize

import (
	"math"

	lbfgsb "github.com/afbarnard/go-lbfgsb"
)

// boundaryDelta is a distance between actual boundary and the
// boundary passed to lbfgsb. This is needed to prevent lbfgs from
// computing likelihood on the boundary, which is sometimes returning
// NAN or infinity.
const boundaryDelta = 1e-5

// LBFGSB is wrapper around go-lbfgsb library. It uses L-BFGS-B
// algorithm of hill climbing.
type LBFGSB struct {
	BaseOptimizer
	dH         float64
	grad       []float64
	exitStatus lbfgsb.ExitStatus
}

// NewLBFGSB creates a new LBFGSB optimizer.
func NewLBFGSB() (lbfgsb *LBFGSB) {
	lbfgsb = &LBFGSB{
		BaseOptimizer: BaseOptimizer{
			repPeriod: 10,
		},
		dH: 1e-6,
	}
	return
}

// Logger is a function which is called back on every iteration.
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

// EvaluateFunction evaluates likelihood function for point x.
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

// EvaluateGradient evaluates gradient for point x.
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
	log.Debugf("Grad(%v)=%v", l.parameters.Values(nil), grad)
	return
}

// Run starts the optimization.
func (l *LBFGSB) Run(iterations int) {
	l.maxL = math.Inf(-1)

	l.PrintHeader()

	bounds := make([][2]float64, len(l.parameters))

	for i, par := range l.parameters {
		bounds[i][0] = par.GetMin() + boundaryDelta
		bounds[i][1] = par.GetMax() - boundaryDelta
	}

	opt := new(lbfgsb.Lbfgsb)
	opt.SetApproximationSize(10)
	opt.SetFTolerance(1e-9)
	opt.SetGTolerance(1e-9)

	opt.SetBounds(bounds)
	opt.SetLogger(l.Logger)

	_, l.exitStatus = opt.Minimize(l, l.parameters.Values(nil))

	switch l.exitStatus.Code {
	case lbfgsb.SUCCESS:
		fallthrough
	case lbfgsb.APPROXIMATE:
		log.Info("LBFGSB status:", l.exitStatus)
	case lbfgsb.WARNING:
		log.Warning("Warning, LBFGSB status:", l.exitStatus)
	default:
		log.Error("Error during LBFGSB:", l.exitStatus)
	}
}

func (l *LBFGSB) Summary() interface{} {
	s := l.BaseOptimizer.Summary().(baseOptimizerSummary)
	s.Status = struct {
		Code       lbfgsb.ExitStatusCode `json:"code"`
		CodeString string                `json:"codeString"`
		Message    string                `json:"message"`
	}{
		l.exitStatus.Code,
		l.exitStatus.Code.String(),
		l.exitStatus.Message,
	}
	return s

}

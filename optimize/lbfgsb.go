package optimize

import (
	"math"

	lbfgsb "github.com/idavydov/go-lbfgsb"
)

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
	l.i++
	err := l.parameters.SetValues(info.X)
	if err != nil {
		panic(err)
	}
	l.PrintLine(l.parameters, -info.F, l.repPeriod)

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

	err := l.parameters.SetValues(x)
	if err != nil {
		panic(err)
	}

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
	err := l.parameters.SetValues(x)
	if err != nil {
		panic(err)
	}

	l1 := -l.Likelihood()
	l.calls++
	for i := range x {
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
	l.SaveStart()

	l.PrintHeader()

	bounds := make([][2]float64, len(l.parameters))

	for i, par := range l.parameters {
		bounds[i][0] = par.GetMin()
		bounds[i][1] = par.GetMax()
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

	l.SaveCheckpoint(true)
	l.saveDeltaT()
}

// Summary returns optimization summary (i.e. success/error, etc).
func (l *LBFGSB) Summary() Summary {
	s := l.BaseOptimizer.Summary().(*baseSummary)
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

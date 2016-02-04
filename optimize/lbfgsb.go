package optimize

import (
	"log"
	"math"

	lbfgsb "github.com/afbarnard/go-lbfgsb"
)

type LBFGSB struct {
	BaseOptimizer
	parameters FloatParameters
	Optimizable
	dH    float64
	grad  []float64
	calls int // likelihood calls
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

func (l *LBFGSB) SetOptimizable(opt Optimizable) {
	l.Optimizable = opt
	l.parameters = opt.GetFloatParameters()
}

func (l *LBFGSB) Logger(info *lbfgsb.OptimizationIterationInformation) {
	l.i = info.Iteration
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
	l.calls += 1
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
	for i, _ := range x {
		no1 := l.Optimizable.Copy()
		par1 := no1.GetFloatParameters()
		par1.SetValues(x)
		v := x[i] - l.dH
		par1[i].Set(v)
		l1 := -no1.Likelihood()
		l.calls += 1

		no2 := no1.Copy()
		par2 := no2.GetFloatParameters()
		v = x[i] + l.dH
		par2[i].Set(v)
		l2 := -no2.Likelihood()
		l.calls += 1

		grad[i] = (l2 - l1) / 2 / l.dH
	}
	select {
	case s := <-l.sig:
		log.Fatal("Received signal exiting:", s)
	default:
	}
	return
}

func (l *LBFGSB) Run(iterations int) {
	l.maxL = math.Inf(-1)
	l.PrintHeader(l.parameters)
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

	log.Print("Exit status: ", exitStatus)

	if !l.Quiet {
		log.Print("Finished LBFGSB")
		log.Printf("Maximum likelihood: %v", l.maxL)
		log.Printf("Likelihood function calls: %v", l.calls)
		log.Printf("Parameter  names: %v", l.parameters.NamesString())
		log.Printf("Parameter values: %v", l.GetMaxLParameters())
	}
	l.PrintFinal(l.parameters)
}

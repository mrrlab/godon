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
	dH   float64
	grad []float64
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
		log.Fatal("Received signal %v, exiting.", s)
	default:
	}
}

func (l *LBFGSB) EvaluateFunction(x []float64) float64 {
	if !l.parameters.ValuesInRange(x) {
		return math.Inf(+1)
	}

	l.parameters.SetValues(x)

	L := l.Likelihood()
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
	if !l.parameters.ValuesInRange(x) {
		for i, par := range l.parameters {
			switch {
			case par.ValueInRange(x[i]):
				grad[i] = 0
			case x[i] < par.GetMin():
				grad[i] = -math.Inf(-1)
			case x[i] > par.GetMax():
				grad[i] = math.Inf(+1)
			default:
				panic("Unknown parameter value")
			}
		}
		return
	}
	no1 := l.Optimizable.Copy()
	par1 := no1.GetFloatParameters()
	par1.SetValues(x)
	l1 := -no1.Likelihood()
	for i, _ := range x {
		no2 := no1.Copy()
		par2 := no2.GetFloatParameters()
		v := x[i] + l.dH
		switch {
		case par2[i].ValueInRange(v):
			par2[i].Set(v)
			l2 := -no2.Likelihood()
			grad[i] = (l2 - l1) / l.dH
		case v < par2[i].GetMin():
			grad[i] = -math.Inf(-1)
		case v > par2[i].GetMax():
			grad[i] = -math.Inf(1)
		default:
			panic("Unknown parameter value")
		}
	}
	return
}

func (l *LBFGSB) Run(iterations int) {
	l.maxL = math.Inf(-1)
	l.PrintHeader(l.parameters)
	bounds := make([][2]float64, len(l.parameters))

	for i, par := range l.parameters {
		//bounds[i] = make([2]float64)
		bounds[i][0] = par.GetMin()
		bounds[i][1] = par.GetMax()
	}

	opt := new(lbfgsb.Lbfgsb)

	opt.SetBounds(bounds)
	opt.SetLogger(l.Logger)

	_, exitStatus := opt.Minimize(l, l.parameters.Values(nil))

	log.Print("Exit status: ", exitStatus)

	if !l.Quiet {
		log.Print("Finished LBFGSB")
		log.Printf("Maximum likelihood: %v", l.maxL)
		log.Printf("Parameter  names: %v", l.parameters.NamesString())
		log.Printf("Parameter values: %v", l.GetMaxLParameters())
	}
	l.PrintFinal(l.parameters)
}

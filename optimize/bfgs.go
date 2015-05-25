package optimize

import (
	"errors"
	"log"
	"math"

	opt "github.com/gonum/optimize"
)

type BFGS struct {
	BaseOptimizer
	parameters FloatParameters
	Optimizable
	dH float64
}

func NewBFGS() (bfgs *BFGS) {
	bfgs = &BFGS{
		BaseOptimizer: BaseOptimizer{
			repPeriod: 10,
		},
		dH: 1e-6,
	}
	return
}

func (b *BFGS) SetOptimizable(opt Optimizable) {
	b.Optimizable = opt
	b.parameters = opt.GetFloatParameters()
}

func (b *BFGS) Init(*opt.FunctionInfo) error {
	return nil
}

func (b *BFGS) Record(l *opt.Location, et opt.EvaluationType, it opt.IterationType, s *opt.Stats) error {
	select {
	case s := <-b.sig:
		log.Printf("Received signal %v, exiting.", s)
		return errors.New("Exiting by signal")
	default:
	}
	return nil
}

func (b *BFGS) Func(x []float64) float64 {
	b.parameters.SetValues(x)
	b.PrintLine(b.parameters, b.Likelihood())
	if !b.parameters.ValuesInRange(x) {
		return math.Inf(+1)
	}

	b.parameters.SetValues(x)

	l := b.Likelihood()
	if l > b.maxL {
		b.maxL = l
		b.maxLPar = b.parameters.Values(b.maxLPar)
	}
	return -l
}

func (b *BFGS) Grad(x, grad []float64) {
	b.parameters.SetValues(x)
	b.PrintLine(b.parameters, b.Likelihood())
	if !b.parameters.ValuesInRange(x) {
		for i, _ := range grad {
			grad[i] = 0
		}
		return
	}
	no1 := b.Optimizable.Copy()
	par1 := no1.GetFloatParameters()
	par1.SetValues(x)
	l1 := -no1.Likelihood()
	for i, _ := range x {
		no2 := no1.Copy()
		par2 := no2.GetFloatParameters()
		par2[i].Set(x[i] + b.dH)
		var l2 float64
		if !par2.InRange() {
			l2 = math.Inf(+1)
		} else {
			l2 = -no2.Likelihood()
		}
		grad[i] = (l2 - l1) / b.dH
	}
	return
}

func (b *BFGS) Run(iterations int) {
	b.maxL = math.Inf(-1)
	b.PrintHeader(b.parameters)
	settings := opt.DefaultSettings()
	settings.MajorIterations = iterations
	settings.Recorder = b

	_, e := opt.Local(b, b.parameters.Values(nil), settings, &opt.BFGS{})
	log.Print(e)

	if !b.Quiet {
		log.Print("Finished BFGS")
		log.Printf("Maximum likelihood: %v", b.maxL)
		log.Printf("Parameter  names: %v", b.parameters.NamesString())
		log.Printf("Parameter values: %v", b.GetMaxLParameters())
	}
	b.PrintFinal(b.parameters)
}

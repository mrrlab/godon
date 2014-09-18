package optimize

import (
	"fmt"
	"log"
	"os"
	"os/signal"
)

type Optimizable interface {
	GetFloatParameters() FloatParameters
	Copy() Optimizable
	Likelihood() float64
}

type Optimizer interface {
	SetOptimizable(Optimizable)
	WatchSignals(...os.Signal)
	SetReportPeriod(period int)
	Run(iterations int)
	GetL() float64
	GetMaxL() float64
	GetMaxLParameters() string
}

type BaseOptimizer struct {
	i         int
	l         float64
	maxL      float64
	maxLPar   []float64
	repPeriod int
	sig       chan os.Signal
	Quiet     bool
}

func (o *BaseOptimizer) WatchSignals(sigs ...os.Signal) {
	o.sig = make(chan os.Signal, 1)
	signal.Notify(o.sig, sigs...)
}

func (o *BaseOptimizer) SetReportPeriod(period int) {
	o.repPeriod = period
}

func (o *BaseOptimizer) PrintHeader(par FloatParameters) {
	if !o.Quiet {
		fmt.Printf("iteration\tlikelihood\t%s\n", par.NamesString())
	}
}

func (o *BaseOptimizer) PrintLine(par FloatParameters, l float64) {
	if !o.Quiet {
		fmt.Printf("%d\t%f\t%s\n", o.i, l, par.ValuesString())
	}
}

func (o *BaseOptimizer) PrintFinal(parameters FloatParameters) {
	if !o.Quiet {
		for _, par := range parameters {
			log.Printf("%s=%v", par.Name(), par.Get())
		}
	}
}

func (o *BaseOptimizer) GetL() float64 {
	return o.l
}
func (o *BaseOptimizer) GetMaxL() float64 {
	return o.maxL
}
func (o *BaseOptimizer) GetMaxLParameters() (s string) {
	for i, v := range o.maxLPar {
		s += fmt.Sprintf("%0.5f", v)
		if i != len(o.maxLPar)-1 {
			s += "\t"
		}
	}
	return s
}

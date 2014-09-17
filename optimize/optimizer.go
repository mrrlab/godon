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
	Optimizable
	parameters FloatParameters
	i          int
	l          float64
	maxL       float64
	maxLPar    string
	repPeriod  int
	sig        chan os.Signal
	Quiet      bool
}

func (o *BaseOptimizer) SetOptimizable(opt Optimizable) {
	o.Optimizable = opt
	o.parameters = opt.GetFloatParameters()
}

func (o *BaseOptimizer) WatchSignals(sigs ...os.Signal) {
	o.sig = make(chan os.Signal, 1)
	signal.Notify(o.sig, sigs...)
}

func (o *BaseOptimizer) SetReportPeriod(period int) {
	o.repPeriod = period
}

func (o *BaseOptimizer) PrintHeader() {
	if !o.Quiet {
		fmt.Printf("iteration\tlikelihood\t%s\n", o.parameters.NamesString())
	}
}

func (o *BaseOptimizer) PrintLine(par FloatParameters, l float64) {
	if !o.Quiet {
		fmt.Printf("%d\t%f\t%s\n", o.i, l, par.ValuesString())
	}
}

func (o *BaseOptimizer) PrintFinal() {
	if !o.Quiet {
		for _, par := range o.parameters {
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
func (o *BaseOptimizer) GetMaxLParameters() string {
	return o.maxLPar
}

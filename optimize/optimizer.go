// Optimize is a collection of optimizers and MCMC samplers.
package optimize

import (
	"fmt"
	"io"
	"os"
	"os/signal"

	"github.com/op/go-logging"
)

// log is the global logging variable.
var log = logging.MustGetLogger("optimize")

// Optimizable is something which can be optimized using the
// optimizer.
type Optimizable interface {
	// GetFloatParameters returns array of float parameters.
	GetFloatParameters() FloatParameters
	// Copy creates a copy of optimizable.
	Copy() Optimizable
	// Likelihood returns likelihood.
	Likelihood() float64
}

// Optimizer an optimizer interface.
type Optimizer interface {
	// SetOptimizable sets model for the optimization.
	SetOptimizable(Optimizable)
	// GetOptimizable returns the underlying model.
	GetOptimizable() Optimizable
	// WatchSignals installs OS hooks to react to signals.
	WatchSignals(...os.Signal)
	// SetReportPeriod specifies how often report line should be printed.
	SetReportPeriod(period int)
	// SetOutput specifies an output writer for the report.
	SetOutput(io.Writer)
	// Starts the optimization or sampling.
	Run(iterations int)
	// GetL return the current likelihood.
	GetL() float64
	// GetMaxL returns the maximum likelihood value.
	GetMaxL() float64
	// GetMaxLParameters returns parameter values for the maximum
	// likelihood value.
	GetMaxLParameters() string
	// GetNCalls return total number of likelihood function calls.
	GetNCalls() int
	// GetNIter returns total number of iterations. Can be smaller than
	// number of likelihood function calls for the gradient-based methods.
	GetNIter() int
	// LoadFromOptimizer can load a model from another optimizer.
	LoadFromOptimizer(Optimizer)
	// PrintFinal prints the final output.
	PrintFinal()
}

// BaseOptimizer contains basic data for an optimizer.
type BaseOptimizer struct {
	Optimizable
	i          int
	calls      int
	l          float64
	maxL       float64
	maxLPar    []float64
	repPeriod  int
	sig        chan os.Signal
	output     io.Writer
	// Quiet controls whether output should be printed.
	Quiet      bool
	parameters FloatParameters
}

// SetOptimizable sets a model for the optimization.
func (o *BaseOptimizer) SetOptimizable(opt Optimizable) {
	o.Optimizable = opt
	o.parameters = opt.GetFloatParameters()
}

// GetOptimizable get the model.
func (o *BaseOptimizer) GetOptimizable() Optimizable {
	return o.Optimizable
}

// LoadFromOptimizer can load a model from another optimizer.
func (o *BaseOptimizer) LoadFromOptimizer(opt Optimizer) {
	o.i = opt.GetNIter()
	o.calls = opt.GetNCalls()
	o.SetOptimizable(opt.GetOptimizable())
}

// WatchSignals installs OS hooks to react to signals.
func (o *BaseOptimizer) WatchSignals(sigs ...os.Signal) {
	o.sig = make(chan os.Signal, 1)
	signal.Notify(o.sig, sigs...)
}

// SetReportPeriod specifies how often report line should be printed.
func (o *BaseOptimizer) SetReportPeriod(period int) {
	o.repPeriod = period
}

// SetOutput specifies an output writer for the report.
func (o *BaseOptimizer) SetOutput(output io.Writer) {
	o.output = output
}

// PrintHeader prints the report header.
func (o *BaseOptimizer) PrintHeader() {
	if !o.Quiet {
		if o.output == nil {
			o.output = os.Stdout
		}
		par := o.Optimizable.GetFloatParameters()
		fmt.Fprintf(o.output, "iteration\tlikelihood\t%s\n", par.NamesString())
	}
}

// PrintLine prints one line of report (iteration, parameter values,
// likelihood).
func (o *BaseOptimizer) PrintLine(par FloatParameters, l float64) {
	if !o.Quiet {
		if o.output == nil {
			o.output = os.Stdout
		}
		fmt.Fprintf(o.output, "%d\t%f\t%s\n", o.i, l, par.ValuesString())
	}
}

// GetNCalls return total number of likelihood function calls.
func (o *BaseOptimizer) GetNCalls() int {
	return o.calls
}

// GetNIter returns total number of iterations. Can be smaller than
// number of likelihood function calls for the gradient-based methods.
func (o *BaseOptimizer) GetNIter() int {
	return o.i
}

// PrintFinal prints the final output.
func (o *BaseOptimizer) PrintFinal() {
	if !o.Quiet {
		par := o.Optimizable.GetFloatParameters()
		log.Noticef("Maximum likelihood: %v", o.maxL)
		log.Infof("Likelihood function calls: %v", o.calls)
		log.Infof("Parameter  names: %v", par.NamesString())
		log.Infof("Parameter values: %v", o.GetMaxLParameters())
		for _, par := range par {
			log.Noticef("%s=%v", par.Name(), par.Get())
		}
	}
}

// GetL return the current likelihood.
func (o *BaseOptimizer) GetL() float64 {
	return o.l
}

// GetMaxL returns the maximum likelihood value. Can me larger than
// current likelihood for MCMC, simulated annealing, etc.
func (o *BaseOptimizer) GetMaxL() float64 {
	return o.maxL
}

// GetMaxLParameters returns parameter values for the maximum
// likelihood value.
func (o *BaseOptimizer) GetMaxLParameters() (s string) {
	for i, v := range o.maxLPar {
		s += fmt.Sprintf("%0.5f", v)
		if i != len(o.maxLPar)-1 {
			s += "\t"
		}
	}
	return s
}

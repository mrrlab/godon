// Package optimize is a collection of optimizers and MCMC samplers.
package optimize

import (
	"fmt"
	"io"
	"os"
	"os/signal"
	"time"

	"github.com/op/go-logging"
)

// log is the global logging variable.
var log = logging.MustGetLogger("optimize")

// Optimizable is something which can be optimized using the
// optimizer.
type Optimizable interface {
	// GetFloatParameters returns array (slice) of float
	// parameters.
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
	// SetTrajectoryOutput specifies an output writer for the trajectory.
	SetTrajectoryOutput(io.Writer)
	// Starts the optimization or sampling.
	Run(iterations int)
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
	// PrintResults prints the results of optimization.
	PrintResults()
	// Summary returns optimization summary for JSON output.
	Summary() Summary
}

// baseSummary stores summary information.
type baseSummary struct {
	// MaxLnL is the maximum likelihood value.
	MaxLnL float64 `json:"maxLnL"`
	// MaxLParameters is the maximum likelihood parameter values.
	MaxLParameters map[string]float64 `json:"maxLParameters"`
	// StartingLnL is the starting likelihood value.
	StartingLnL float64 `json:"startingLnL"`
	// StartingParameters is the starting parameter values.
	StartingParameters map[string]float64 `json:"startingParameters"`
	// NIterations is the number of iterations, for gradient-based methods can be less than number of funciton calls.
	NIterations int `json:"nIterations"`
	// NCalls is the number of function calls, cumulative for chained methods
	NCalls int `json:"nLikelihoodComputations"`
	// Status is the optimization status.
	Status interface{} `json:"status,omitempty"`
	// OptimizationTime is the optimization time in seconds.
	OptimizationTime float64 `json:"optimizationTime,omitempty"`
}

// Summary allows quering of maximum likelihood estimates.
type Summary interface {
	GetMaxLikelihood() float64
	GetMaxLikelihoodParameters() map[string]float64
}

func (b *baseSummary) GetMaxLikelihood() float64 {
	return b.MaxLnL
}

func (b *baseSummary) GetMaxLikelihoodParameters() map[string]float64 {
	r := make(map[string]float64, len(b.MaxLParameters))
	for k, v := range b.MaxLParameters {
		r[k] = v
	}
	return r
}

// BaseOptimizer contains basic data for an optimizer.
type BaseOptimizer struct {
	Optimizable
	i         int
	calls     int
	maxL      float64
	maxLPar   []float64
	startL    float64
	startPar  []float64
	repPeriod int
	sig       chan os.Signal
	output    io.Writer
	// Quiet controls whether output should be printed.
	Quiet      bool
	parameters FloatParameters

	otime, ftime float64
	//startTime is a starting time for deltaT compuataions.
	startTime time.Time
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

// SetTrajectoryOutput specifies an output writer for the trajectory.
func (o *BaseOptimizer) SetTrajectoryOutput(output io.Writer) {
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

// PrintResults prints the optimization results.
func (o *BaseOptimizer) PrintResults() {
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

// GetMaxL returns the maximum likelihood value. Can me larger than
// current likelihood for MCMC, simulated annealing, etc.
func (o *BaseOptimizer) GetMaxL() float64 {
	return o.maxL
}

// GetMaxLParameters returns parameter values for the maximum
// likelihood value as a string.
func (o *BaseOptimizer) GetMaxLParameters() (s string) {
	for i, v := range o.maxLPar {
		s += fmt.Sprintf("%0.5f", v)
		if i != len(o.maxLPar)-1 {
			s += "\t"
		}
	}
	return s
}

// GetParametersMap returns parameter values as a map.
func (o *BaseOptimizer) GetParametersMap(par []float64) (m map[string]float64) {
	m = make(map[string]float64, len(par))
	for i, v := range par {
		m[o.parameters[i].Name()] = v
	}
	return m
}

// saveDeltaT saves computations time in seconds.
func (o *BaseOptimizer) saveDeltaT() {
	o.otime = deltaT(o.startTime)
}

// SaveStart saves starting point and likelihood before optimization.
func (o *BaseOptimizer) SaveStart() {
	o.startTime = time.Now()
	l := o.Likelihood()
	o.calls++
	o.startL = l
	o.startPar = o.parameters.Values(nil)
	o.maxL = l
	o.maxLPar = o.parameters.Values(nil)
}

// Summary returns optimization summary.
func (o *BaseOptimizer) Summary() Summary {
	return &baseSummary{
		MaxLnL:             o.maxL,
		MaxLParameters:     o.GetParametersMap(o.maxLPar),
		StartingLnL:        o.startL,
		StartingParameters: o.GetParametersMap(o.startPar),
		NIterations:        o.GetNIter(),
		NCalls:             o.GetNCalls(),
		OptimizationTime:   o.otime,
	}
}

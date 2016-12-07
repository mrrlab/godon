package main

import (
	"fmt"
	"os"

	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/optimize"
)

// optimizerSettings stores settings for creation of a new optimizer.
type optimizerSettings struct {
	method string
	model  cmodel.TreeOptimizableSiteClass

	iterations int

	report int

	accept   int
	adaptive bool
	skip     int
	maxAdapt int

	trajF *os.File

	seed int64
}

// newOptimzerSettings creates a new optimizerSettings from
// the command line parameters (global variables).
func newOptimzerSettings(model cmodel.TreeOptimizableSiteClass) *optimizerSettings {
	return &optimizerSettings{
		method: *method,
		model:  model,

		iterations: *iterations,

		report: *report,

		accept:   *accept,
		adaptive: *adaptive,
		skip:     *skip,
		maxAdapt: *maxAdapt,

		trajF: trajF,

		seed: *seed,
	}
}

// create creates and initializes a new optimizer from optimizerSettings.
func (o *optimizerSettings) create() (optimize.Optimizer, error) {
	// iteration to skip before annealing, for adaptive mcmc
	if o.adaptive {
		as := optimize.NewAdaptiveSettings()
		if o.skip < 0 {
			o.skip = o.iterations / 20
		}
		if o.maxAdapt < 0 {
			o.maxAdapt = o.iterations / 5
		}
		log.Infof("Setting adaptive parameters, skip=%v, maxAdapt=%v", o.skip, o.maxAdapt)
		as.Skip = o.skip
		as.MaxAdapt = o.maxAdapt
		o.model.SetAdaptive(as)
	}

	opt, err := o.getOptimizer()
	if err != nil {
		return nil, err
	}
	log.Infof("Using %s optimization.", o.method)

	opt.SetTrajectoryOutput(o.trajF)
	opt.SetOptimizable(o.model)

	opt.SetReportPeriod(o.report)

	return opt, nil
}

// getOptimizer returns an optimizer from settings.
func (o *optimizerSettings) getOptimizer() (optimize.Optimizer, error) {
	switch o.method {
	case "lbfgsb":
		return optimize.NewLBFGSB(), nil
	case "simplex":
		return optimize.NewDS(), nil
	case "mh":
		chain := optimize.NewMH(false, 0)
		chain.AccPeriod = o.accept
		return chain, nil
	case "annealing":
		chain := optimize.NewMH(true, o.maxAdapt)
		chain.AccPeriod = o.accept
		return chain, nil
	case "n_lbfgs":
		return optimize.NewNLOPT(optimize.NLOPT_LBFGS, o.seed), nil
	case "n_simplex":
		return optimize.NewNLOPT(optimize.NLOPT_SIMPLEX, o.seed), nil
	case "n_cobyla":
		return optimize.NewNLOPT(optimize.NLOPT_COBYLA, o.seed), nil
	case "n_bobyqa":
		return optimize.NewNLOPT(optimize.NLOPT_BOBYQA, o.seed), nil
	case "n_sqp":
		return optimize.NewNLOPT(optimize.NLOPT_SQP, o.seed), nil
	case "n_direct":
		return optimize.NewNLOPT(optimize.NLOPT_DIRECT, o.seed), nil
	case "n_crs":
		return optimize.NewNLOPT(optimize.NLOPT_CRS, o.seed), nil
	case "n_mlsl":
		return optimize.NewNLOPT(optimize.NLOPT_MLSL, o.seed), nil
	case "none":
		return optimize.NewNone(), nil
	}
	return nil, fmt.Errorf("Unknown optimization method: %s", o.method)
}

package main

import (
	"fmt"
	"os"
	"time"

	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/optimize"
)

// getOptimizerFromString returns an optimizer from a string.
func getOptimizerFromString(method string, accept, annealingSkip int, seed int64) (optimize.Optimizer, error) {
	switch method {
	case "lbfgsb":
		return optimize.NewLBFGSB(), nil
	case "simplex":
		return optimize.NewDS(), nil
	case "mh":
		chain := optimize.NewMH(false, 0)
		chain.AccPeriod = accept
		return chain, nil
	case "annealing":
		chain := optimize.NewMH(true, annealingSkip)
		chain.AccPeriod = accept
		return chain, nil
	case "n_lbfgs":
		return optimize.NewNLOPT(optimize.NLOPT_LBFGS, seed), nil
	case "n_simplex":
		return optimize.NewNLOPT(optimize.NLOPT_SIMPLEX, seed), nil
	case "n_cobyla":
		return optimize.NewNLOPT(optimize.NLOPT_COBYLA, seed), nil
	case "n_bobyqa":
		return optimize.NewNLOPT(optimize.NLOPT_BOBYQA, seed), nil
	case "n_sqp":
		return optimize.NewNLOPT(optimize.NLOPT_SQP, seed), nil
	case "n_direct":
		return optimize.NewNLOPT(optimize.NLOPT_DIRECT, seed), nil
	case "n_crs":
		return optimize.NewNLOPT(optimize.NLOPT_CRS, seed), nil
	case "n_mlsl":
		return optimize.NewNLOPT(optimize.NLOPT_MLSL, seed), nil
	case "none":
		return optimize.NewNone(), nil
	}
	return nil, fmt.Errorf("Unknown optimization method: %s", method)
}

func newData() (*cmodel.Data, error) {
	data, err := cmodel.NewData(*gcodeID, *alignmentFileName, *treeFileName, *cFreq)

	if err != nil {
		return nil, err
	}

	if len(*cFreqFileName) > 0 {
		data.SetCodonFreqFromFile(*cFreqFileName)
	}

	if *fgBranch >= 0 {
		data.SetForegroundBranch(*fgBranch)
	}

	if *model == "BS" || *model == "BSG" || *model == "BSGE" {
		if data.GetNClass1() < 1 {
			log.Warning("Warning: no class=1 nodes")
		}
	}

	return data, nil
}

func runOptimization(h0 bool, start map[string]float64) (summary OptimizationSummary) {
	startTime := time.Now()

	data, err := newData()
	if err != nil {
		log.Fatal(err)
	}

	summary.StartingTree = data.Tree.ClassString()

	ms := newModelSettings(data)

	m, err := ms.createInitalized()
	if err != nil {
		log.Fatal(err)
	}

	if len(start) > 0 {
		par := m.GetFloatParameters()
		err := par.SetFromMap(start)
		if err != nil {
			log.Fatal(err)
		}
	}

	// iteration to skip before annealing, for adaptive mcmc
	annealingSkip := 0
	if *adaptive {
		as := optimize.NewAdaptiveSettings()
		if *skip < 0 {
			*skip = *iterations / 20
		}
		if *maxAdapt < 0 {
			*maxAdapt = *iterations / 5
		}
		annealingSkip = *maxAdapt
		log.Infof("Setting adaptive parameters, skip=%v, maxAdapt=%v", *skip, *maxAdapt)
		as.Skip = *skip
		as.MaxAdapt = *maxAdapt
		m.SetAdaptive(as)
	}

	f := os.Stdout

	if *trajF != "" {
		f, err = os.OpenFile(*trajF, os.O_WRONLY|os.O_CREATE|os.O_APPEND, 0666)
		if err != nil {
			log.Fatal("Error creating trajectory file:", err)
		}
		defer f.Close()
	}

	var opt optimize.Optimizer

	opt, err = getOptimizerFromString(*method, *accept, annealingSkip, *seed)
	if err != nil {
		log.Fatal(err)
	}

	log.Infof("Using %s optimization.", *method)

	opt.SetTrajectoryOutput(f)
	opt.SetOptimizable(m)

	opt.SetReportPeriod(*report)

	opt.Run(*iterations)
	summary.Optimizer = opt.Summary()

	opt.PrintResults()
	if !*noFinal {
		m.Final()
	}
	summary.Model = m.Summary()

	if !*noOptBrLen {
		log.Infof("outtree=%s", data.Root())
		data.Root()
		summary.FinalTree = data.Tree.ClassString()
	}

	if *outTreeF != "" {
		f, err := os.Create(*outTreeF)
		if err != nil {
			log.Error("Error creating tree output file:", err)
		} else {
			f.WriteString(data.Tree.String() + "\n")
			f.Close()
		}
	}

	if *printFull && (ms.aggMode != cmodel.AggNone) {
		m.SetAggregationMode(cmodel.AggNone)
		L := m.Likelihood()
		log.Notice("Full likelihood: ", L)
		summary.FullLnL = L
	}

	endTime := time.Now()

	deltaT := endTime.Sub(startTime)
	summary.Time = deltaT.Seconds()

	return
}

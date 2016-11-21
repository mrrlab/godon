package main

import (
	"errors"
	"fmt"
	"os"
	"time"

	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/optimize"
)

// getAggModeFromString returns an aggregation mode constant from cmodel
// from a string.
func getAggModeFromString(aggModeString string) (cmodel.AggMode, error) {
	switch aggModeString {
	case "none":
		return cmodel.AggNone, nil
	case "observed":
		return cmodel.AggObserved, nil
	case "observed_new":
		return cmodel.AggObservedNew, nil
	case "fixed":
		return cmodel.AggFixed, nil
	case "random":
		return cmodel.AggRandom, nil
	}
	return cmodel.AggNone, fmt.Errorf("Unknown aggregation mode: %s", aggModeString)
}

// getModelFromString returns a model from string and other parameters.
func getModelFromString(model string, data *cmodel.Data, fixw bool,
	ncatb, ncatsg, ncatcg int) (cmodel.TreeOptimizableSiteClass, error) {
	switch model {
	case "M0":
		log.Info("Using M0 model")
		return cmodel.NewM0(data), nil
	case "M1a":
		log.Info("Using M1a model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewM2(data, false, ncatsg, ncatcg), nil
	case "M2a":
		log.Info("Using M2a model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewM2(data, true, ncatsg, ncatcg), nil
	case "M7":
		log.Info("Using M7 model")
		log.Infof("%d beta categories, %d site gamma categories, %d codon gama categories", ncatb, ncatsg, ncatcg)
		return cmodel.NewM8(data, false, false, ncatb, ncatsg, ncatcg), nil
	case "M8":
		log.Info("Using M8 model")
		log.Infof("%d beta categories, %d site gamma categories, %d codon gama categories", ncatb, ncatsg, ncatcg)
		return cmodel.NewM8(data, true, fixw, ncatb, ncatsg, ncatcg), nil
	case "BSG":
		log.Info("Using branch site gamma model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewBranchSiteGamma(data, fixw, ncatsg, ncatcg), nil
	case "BSGE":
		log.Info("Using branch site gamma model with explicit rates")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewBranchSiteGammaERates(data, fixw, ncatsg, ncatcg), nil
	case "BS":
		log.Info("Using branch site model")
		return cmodel.NewBranchSite(data, fixw), nil
	}
	return nil, errors.New("Unknown model specification")
}

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

func runOptimization(h0 bool, start map[string]float64) (summary OptimizationSummary) {
	startTime := time.Now()

	data, err := cmodel.NewData(*gcodeID, *alignmentFileName, *treeFileName, *cFreq)

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

	summary.StartingTree = data.Tree.ClassString()
	if err != nil {
		log.Fatal(err)
	}

	m, err := getModelFromString(*model, data, h0, *ncatb, *ncatsg, *ncatcg)
	if err != nil {
		log.Fatal(err)
	}

	log.Infof("Model has %d site class(es)", m.GetNClass())

	if !*noOptBrLen {
		log.Info("Will optimize branch lengths")
		log.Infof("Maximum branch length: %f", *maxBrLen)
		m.SetMaxBranchLength(*maxBrLen)
		m.SetOptimizeBranchLengths()
	} else {
		log.Info("Will not optimize branch lengths")
	}

	aggMode, err := getAggModeFromString(*aggregate)
	if err != nil {
		log.Fatal(err)
	}
	if aggMode != cmodel.AggNone {
		log.Infof("Aggregation mode: %s", *aggregate)
	}
	m.SetAggregationMode(aggMode)

	if len(start) > 0 {
		par := m.GetFloatParameters()
		err := par.SetFromMap(start)
		if err != nil {
			log.Fatal(err)
		}
	} else if *startF != "" {
		l, err := lastLine(*startF)
		par := m.GetFloatParameters()
		if err == nil {
			err = par.ReadLine(l)
		}
		if err != nil {
			log.Debug("Reading start file as JSON")
			err2 := par.ReadFromJSON(*startF)
			// startF is neither trajectory nor correct JSON
			if err2 != nil {
				log.Error("Error reading start position from JSON:", err2)
				log.Fatal("Error reading start position from trajectory file:", err)
			}
		}
		if !par.InRange() {
			log.Fatal("Initial parameters are not in the range")
		}
	} else if *randomize {
		log.Info("Using uniform (in the boundaries) random starting point")
		par := m.GetFloatParameters()
		par.Randomize()
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

	log.Infof("Model has %d parameters.", len(m.GetFloatParameters()))

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

	if *printFull && (aggMode != cmodel.AggNone) {
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

package main

import (
	"errors"
	"fmt"
	"os"
	"time"

	"bitbucket.org/Davydov/godon/bio"
	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/tree"
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
func getModelFromString(model string, cali codon.Sequences, t *tree.Tree, cf codon.Frequency,
	fixw bool, ncatb, ncatsg, ncatcg int) (cmodel.TreeOptimizableSiteClass, error) {
	switch model {
	case "M0":
		log.Info("Using M0 model")
		return cmodel.NewM0(cali, t, cf), nil
	case "M1a":
		log.Info("Using M1a model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewM2(cali, t, cf, false, ncatsg, ncatcg), nil
	case "M2a":
		log.Info("Using M2a model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewM2(cali, t, cf, true, ncatsg, ncatcg), nil
	case "M7":
		log.Info("Using M7 model")
		log.Infof("%d beta categories, %d site gamma categories, %d codon gama categories", ncatb, ncatsg, ncatcg)
		return cmodel.NewM8(cali, t, cf, false, false, ncatb, ncatsg, ncatcg), nil
	case "M8":
		log.Info("Using M8 model")
		log.Infof("%d beta categories, %d site gamma categories, %d codon gama categories", ncatb, ncatsg, ncatcg)
		return cmodel.NewM8(cali, t, cf, true, fixw, ncatb, ncatsg, ncatcg), nil
	case "BSG":
		log.Info("Using branch site gamma model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewBranchSiteGamma(cali, t, cf, fixw, ncatsg, ncatcg), nil
	case "BSGE":
		log.Info("Using branch site gamma model with explicit rates")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewBranchSiteGammaERates(cali, t, cf, fixw, ncatsg, ncatcg), nil
	case "BS":
		log.Info("Using branch site model")
		return cmodel.NewBranchSite(cali, t, cf, fixw), nil
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

func runOptimization(startFileName string, h0 bool) (summary *RunSummary) {
	startTime := time.Now()
	summary = &RunSummary{}

	gcode, ok := bio.GeneticCodes[*gcodeID]
	if !ok {
		log.Fatalf("couldn't load genetic code with id=%d", gcodeID)
	}
	log.Infof("Genetic code: %d, \"%s\"", gcode.ID, gcode.Name)
	fastaFile, err := os.Open(*alignmentFileName)
	if err != nil {
		log.Fatal(err)
	}
	defer fastaFile.Close()

	ali, err := bio.ParseFasta(fastaFile)
	if err != nil {
		log.Fatal(err)
	}

	cali, err := codon.ToCodonSequences(ali, gcode)
	if err != nil {
		log.Fatal(err)
	}

	if cali.Length() == 0 {
		log.Fatal("Zero length alignment")
	}
	log.Infof("Read alignment of %d codons, %d fixed positions, %d ambiguous positions", cali.Length(), cali.NFixed(), cali.NAmbiguous())

	treeFile, err := os.Open(*treeFileName)
	if err != nil {
		log.Fatal(err)
	}
	defer treeFile.Close()

	t, err := tree.ParseNewick(treeFile)
	if err != nil {
		log.Fatal(err)
	}

	log.Debugf("intree=%s", t)
	log.Debugf("brtree=%s", t.BrString())

	// Root the tree in the end
	var root = false
	var rootID = 0

	if t.IsRooted() {
		log.Warning("Tree is rooted. Will unroot.")
		root = true
		rootID, err = t.Unroot()
		if err != nil {
			log.Fatal("Error unrooting tree:", err)
		}
	}

	log.Infof("intree_unroot=%s", t)
	summary.StartingTree = t.ClassString()
	log.Debugf("brtree_unroot=%s", t.BrString())
	log.Debug(t.FullString())

	var cf codon.Frequency

	if *cFreqFileName != "" {
		cFreqFile, err := os.Open(*cFreqFileName)
		if err != nil {
			log.Fatal(err)
		}
		cf, err = codon.ReadFrequency(cFreqFile, gcode)
		if err != nil {
			log.Fatal(err)
		}
	} else {
		switch *cFreq {
		case "F0":
			log.Info("F0 frequency")
			cf = codon.F0(cali)
		case "F3X4":
			log.Info("F3X4 frequency")
			cf = codon.F3X4(cali)
		default:
			log.Fatal("Unknow codon freuquency specification")
		}
	}
	log.Debug(cf)

	if *fgBranch >= 0 {
		for _, node := range t.NodeIDArray() {
			if node == nil {
				continue
			}
			if node.ID == *fgBranch {
				node.Class = 1
			} else {
				node.Class = 0
			}
		}
	} else {
		class1 := 0
		for range t.ClassNodes(1) {
			class1++
		}
		if class1 == 0 &&
			(*model == "BS" || *model == "BSG" || *model == "BSGE") {
			log.Warning("Warning: no class=1 nodes")
		}
	}

	m, err := getModelFromString(*model, cali, t, cf, h0, *ncatb, *ncatsg, *ncatcg)
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

	if startFileName != "" {
		l, err := lastLine(startFileName)
		par := m.GetFloatParameters()
		if err == nil {
			err = par.ReadLine(l)
		}
		if err != nil {
			log.Debug("Reading start file as JSON")
			err2 := par.ReadFromJSON(startFileName)
			// startFileName is neither trajectory nor correct JSON
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

	if *outF != "" {
		f, err = os.Create(*outF)
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

	opt.SetOutput(f)
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
		if root {
			log.Infof("unrooted_outtree=%s", t)
			err = t.Root(rootID)
			if err != nil {
				log.Error("Error rooting tree:", err)
			}
		}
		log.Infof("outtree=%s", t)
		summary.FinalTree = t.ClassString()
	}

	if *outTreeF != "" {
		f, err := os.Create(*outTreeF)
		if err != nil {
			log.Error("Error creating tree output file:", err)
		} else {
			f.WriteString(t.String() + "\n")
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
	log.Noticef("Running time: %v", deltaT)
	summary.Time = deltaT.Seconds()

	return
}

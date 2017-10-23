package main

import (
	"os"

	"bitbucket.org/Davydov/godon/cmodel"
)

// newData reads data using command line parameters.
func newData() (*cmodel.Data, error) {
	data, err := cmodel.NewData(*gcodeID, *alignmentFileName, *treeFileName, *cFreq)

	if err != nil {
		return nil, err
	}

	if !*noOptBrLen {
		err = data.Unroot()
		if err != nil {
			return nil, err
		}
	}

	if len(*cFreqFileName) > 0 {
		err := data.SetCodonFreqFromFile(*cFreqFileName)
		if err != nil {
			return nil, err
		}
	}

	if *fgBranch >= 0 {
		data.SetForegroundBranch(*fgBranch)
	}

	if *model == "BS" || *model == "BSG" || *model == "BSGE" {
		if data.GetNClass1() < 1 {
			log.Warning("No class=1 nodes (#1)")
		}
	}

	return data, nil
}

// runOptimization runs optimization for model with optimizers
// settings and optional starting point. If minLikelihood > 0,
// optimization is always performed. If 0 >= minLikelihood >
// startingLikelihood, no optimization is performed.
func runOptimization(m cmodel.TreeOptimizableSiteClass, o *optimizerSettings, start map[string]float64, minLikelihood float64, quiet bool) (summary OptimizationSummary) {
	if m.GetOptimizeBranchLengths() {
		summary.StartingTree = m.GetTreeString()
	}

	if len(start) > 0 {
		setStart(m, start)
	}

	if minLikelihood <= 0 && m.Likelihood() < minLikelihood {
		// restore method before leaving the function
		defer func(savedMethod string) {
			o.method = savedMethod
		}(o.method)

		log.Info("Won't optimize, since starting point is no better")
		o.method = "none"
	}

	opt, err := o.create()
	if err != nil {
		log.Fatal(err)
	}

	opt.Run(o.iterations)
	summary.Optimizer = opt.Summary()

	opt.PrintResults(quiet)

	if m.GetOptimizeBranchLengths() {
		summary.FinalTree = m.GetTreeString()
	}

	return summary
}

// setStart sets starting point for model, or logs error message and exits.
func setStart(m cmodel.TreeOptimizableSiteClass, start map[string]float64) {
	par := m.GetFloatParameters()
	err := par.SetFromMap(start)
	if err != nil {
		log.Fatal(err)
	}
}

// computeFinal computes BEB & NEB for a given point
func computeFinal(m cmodel.TreeOptimizableSiteClass, start map[string]float64) interface{} {
	setStart(m, start)
	m.Final(*neb, *beb, *codonRates, *siteRates, *codonOmega)
	return m.Summary()
}

// optimization is calling optimization.
func optimization() OptimizationSummary {
	data, err := newData()
	if err != nil {
		log.Fatal(err)
	}

	ms := newModelSettings(data)

	m, err := ms.createInitalized(false)
	if err != nil {
		log.Fatal(err)
	}

	o := newOptimizerSettings(m)

	summary := runOptimization(m, o, nil, 1, false)

	if *final {
		m.Final(*neb, *beb, *codonRates, *siteRates, *codonOmega)
	}
	summary.Model = m.Summary()

	if !*noOptBrLen {
		err := data.Root()
		if err != nil {
			log.Error(err)
		}
		log.Infof("outtree=%s", data.Tree)
	}

	if *outTreeF != "" {
		f, err := os.Create(*outTreeF)
		if err != nil {
			log.Error("Error creating tree output file:", err)
		} else {
			_, err := f.WriteString(data.Tree.ClassString() + "\n")
			if err != nil {
				log.Error(err)
			}
			f.Close()
		}
	}

	if *printFull && (ms.aggMode != cmodel.AggNone) {
		m.SetAggregationMode(cmodel.AggNone)
		L := m.Likelihood()
		log.Notice("Full likelihood: ", L)
		summary.FullLnL = L
	}

	return summary
}

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
			log.Warning("Warning: no class=1 nodes")
		}
	}

	return data, nil
}

// runOptimization runs optimization for model with optimizers settings
// and optional starting point.
func runOptimization(m cmodel.TreeOptimizableSiteClass, o *optimizerSettings, start map[string]float64) (summary OptimizationSummary) {
	if m.GetOptimizeBranchLengths() {
		summary.StartingTree = m.GetTreeString()
	}

	opt, err := o.create()
	if err != nil {
		log.Fatal(err)
	}

	if len(start) > 0 {
		setStart(m, start)
	}

	opt.Run(o.iterations)
	summary.Optimizer = opt.Summary()

	opt.PrintResults()

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
	m.Final(*neb, *beb, *codonRates, *codonOmega)
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

	summary := runOptimization(m, o, nil)

	if *final {
		m.Final(*neb, *beb, *codonRates, *codonOmega)
	}
	summary.Model = m.Summary()

	if !*noOptBrLen {
		log.Infof("outtree=%s", data.Root())
		err := data.Root()
		if err != nil {
			log.Error(err)
		}
	}

	if *outTreeF != "" {
		f, err := os.Create(*outTreeF)
		if err != nil {
			log.Error("Error creating tree output file:", err)
		} else {
			_, err := f.WriteString(data.Tree.StringPrecise() + "\n")
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

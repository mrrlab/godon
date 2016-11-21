package main

import (
	"os"
	"time"

	"bitbucket.org/Davydov/godon/cmodel"
)

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

	o := newOptimzerSettings(m)
	opt, err := o.create()
	if err != nil {
		log.Fatal(err)
	}

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

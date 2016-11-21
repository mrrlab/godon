package main

import (
	"errors"
	"fmt"

	"bitbucket.org/Davydov/godon/cmodel"
)

type modelSettings struct {
	name   string
	data   *cmodel.Data
	fixw   bool
	ncatb  int
	ncatsg int
	ncatcg int

	noOptBrLen  bool
	maxBrLen    float64
	aggModeName string
	aggMode     cmodel.AggMode

	startF    string
	randomize bool
}

func newModelSettings(data *cmodel.Data) *modelSettings {
	return &modelSettings{
		name:   *model,
		data:   data,
		fixw:   *fixw,
		ncatb:  *ncatb,
		ncatsg: *ncatsg,
		ncatcg: *ncatcg,

		noOptBrLen:  *noOptBrLen,
		maxBrLen:    *maxBrLen,
		aggModeName: *aggregate,

		startF:    *startF,
		randomize: *randomize,
	}
}

func (ms *modelSettings) createModel() (cmodel.TreeOptimizableSiteClass, error) {
	switch ms.name {
	case "M0":
		log.Info("Using M0 model")
		return cmodel.NewM0(ms.data), nil
	case "M1a":
		log.Info("Using M1a model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewM2(ms.data, false, ms.ncatsg, ms.ncatcg), nil
	case "M2a":
		log.Info("Using M2a model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewM2(ms.data, true, ms.ncatsg, ms.ncatcg), nil
	case "M7":
		log.Info("Using M7 model")
		log.Infof("%d beta categories, %d site gamma categories, %d codon gama categories", ncatb, ncatsg, ncatcg)
		return cmodel.NewM8(ms.data, false, false, ms.ncatb, ms.ncatsg, ms.ncatcg), nil
	case "M8":
		log.Info("Using M8 model")
		log.Infof("%d beta categories, %d site gamma categories, %d codon gama categories", ncatb, ncatsg, ncatcg)
		return cmodel.NewM8(ms.data, true, ms.fixw, ms.ncatb, ms.ncatsg, ms.ncatcg), nil
	case "BSG":
		log.Info("Using branch site gamma model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewBranchSiteGamma(ms.data, ms.fixw, ms.ncatsg, ms.ncatcg), nil
	case "BSGE":
		log.Info("Using branch site gamma model with explicit rates")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewBranchSiteGammaERates(ms.data, ms.fixw, ms.ncatsg, ms.ncatcg), nil
	case "BS":
		log.Info("Using branch site model")
		return cmodel.NewBranchSite(ms.data, ms.fixw), nil
	}
	return nil, errors.New("Unknown model specification")
}

// getAggModereturns an aggregation mode constant from cmodel
// from a string.
func (ms *modelSettings) getAggMode() (cmodel.AggMode, error) {
	switch ms.aggModeName {
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
	return cmodel.AggNone, fmt.Errorf("Unknown aggregation mode: %s", ms.aggModeName)
}

func (ms *modelSettings) createInitalized() (cmodel.TreeOptimizableSiteClass, error) {
	m, err := ms.createModel()

	log.Infof("Model has %d parameters.", len(m.GetFloatParameters()))

	if err != nil {
		return nil, err
	}

	if !ms.noOptBrLen {
		log.Info("Will optimize branch lengths")
		log.Infof("Maximum branch length: %f", ms.maxBrLen)
		m.SetMaxBranchLength(ms.maxBrLen)
		m.SetOptimizeBranchLengths()
	} else {
		log.Info("Will not optimize branch lengths")
	}

	ms.aggMode, err = ms.getAggMode()
	if err != nil {
		return nil, err
	}

	if ms.aggMode != cmodel.AggNone {
		log.Infof("Aggregation mode: %s", ms.aggModeName)
	}
	m.SetAggregationMode(ms.aggMode)

	if ms.startF != "" {
		l, err := lastLine(ms.startF)
		par := m.GetFloatParameters()
		if err == nil {
			err = par.ReadLine(l)
		}
		if err != nil {
			log.Debug("Reading start file as JSON")
			err2 := par.ReadFromJSON(ms.startF)
			// startF is neither trajectory nor correct JSON
			if err2 != nil {
				log.Error("Error reading start position from JSON:", err2)
				return nil, fmt.Errorf("Error reading start position from trajectory file: %v", err)
			}
		}
		if !par.InRange() {
			return nil, errors.New("Initial parameters are not in the range")
		}
	} else if ms.randomize {
		log.Info("Using uniform (in the boundaries) random starting point")
		par := m.GetFloatParameters()
		par.Randomize()
	}

	log.Infof("Model has %d site class(es)", m.GetNClass())

	return m, nil
}

package main

import (
	"testing"

	"github.com/op/go-logging"

	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/optimize"
)

func init() {
	// disable logging for benchmarks
	logging.SetLevel(logging.WARNING, "optimize")
	logging.SetLevel(logging.ERROR, "cmodel")
}

func TestBSCG(tst *testing.T) {
	data, err := cmodel.GetTreeAlignment(data1, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	bsg := cmodel.NewBranchSiteGamma(data, true, 1, 4)
	m := optimize.Optimizable(bsg).Copy()
	npar := len(bsg.GetFloatParameters())
	if npar != 5 {
		tst.Error("Wrong number of parameters for BSG:", npar)
	}

	ds := optimize.NewDS()
	ds.SetOptimizable(m)
	ds.Quiet = true
	ds.Run(5)

	npar1 := len(m.GetFloatParameters())
	m = m.Copy()
	npar2 := len(m.GetFloatParameters())
	if npar != npar1 || npar1 != npar2 {
		tst.Error("Parameter number mismatch after copy:", npar, npar1, npar2)
	}

	mh := optimize.NewMH(false, 0)
	mh.SetOptimizable(m)
	mh.Quiet = true
	mh.Run(5)
}

func TestBSSG(tst *testing.T) {
	data, err := cmodel.GetTreeAlignment(data1, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	bsg := cmodel.NewBranchSiteGamma(data, true, 4, 1)
	m := optimize.Optimizable(bsg).Copy()
	npar := len(bsg.GetFloatParameters())
	if npar != 5 {
		tst.Error("Wrong number of parameters for BSG:", npar)
	}

	ds := optimize.NewDS()
	ds.SetOptimizable(m)
	ds.Quiet = true
	ds.Run(5)

	npar1 := len(m.GetFloatParameters())
	m = m.Copy()
	npar2 := len(m.GetFloatParameters())
	if npar != npar1 || npar1 != npar2 {
		tst.Error("Parameter number mismatch after copy:", npar, npar1, npar2)
	}

	mh := optimize.NewMH(false, 0)
	mh.SetOptimizable(m)
	mh.Quiet = true
	mh.Run(5)
}

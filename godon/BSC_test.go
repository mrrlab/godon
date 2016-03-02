package main

import (
	"testing"

	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/optimize"
)

func TestBSC(tst *testing.T) {
	t, cali, err := cmodel.GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := cmodel.F3X4(cali)
	bsc := cmodel.NewBranchSiteC(cali, t, cf)
	m := optimize.Optimizable(bsc).Copy()
	npar := len(bsc.GetFloatParameters())
	if npar != 4 {
		tst.Error("Wrong number of parameters for BSC:", npar)
	}

	ds := optimize.NewDS()
	ds.SetOptimizable(m)
	ds.Run(5)

	npar1 := len(m.GetFloatParameters())
	m = m.Copy()
	npar2 := len(m.GetFloatParameters())
	if npar != npar1 || npar1 != npar2 {
		tst.Error("Parameter number mismatch after copy:", npar, npar1, npar2)
	}

	mh := optimize.NewMH(false, 0)
	mh.SetOptimizable(m)
	mh.Run(5)
}

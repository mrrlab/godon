package main

import (
	"testing"

	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/optimize"
)

func TestBSG(tst *testing.T) {
	t, cali, err := cmodel.GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := cmodel.F3X4(cali)
	bsg := cmodel.NewBranchSiteGamma(cali, t, cf, 4, true)
	m := optimize.Optimizable(bsg).Copy()
	npar := len(bsg.GetFloatParameters())
	if npar != 5 {
		tst.Error("Wrong number of parameters for BSG:", npar)
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

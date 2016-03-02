package main

import (
	"testing"

	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/optimize"
)

const (
	data1 = "EMGT00050000008747.Drosophila.002"
	data2 = "ENSGT00550000073950.Euteleostomi.07.001"
	data3 = "EMGT00050000000025.Drosophila.001"
)

func TestSimplex(tst *testing.T) {
	t, cali, err := cmodel.GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := cmodel.F3X4(cali)
	m0 := cmodel.NewM0(cali, t, cf)
	m0.SetParameters(2, 0.5)
	npar := len(m0.GetFloatParameters())
	if npar != 2 {
		tst.Error("Wrong number of parameters for M0:", npar)
	}

	ds := optimize.NewDS()
	ds.SetOptimizable(m0)
	ds.Run(5)
}

func TestMH(tst *testing.T) {
	t, cali, err := cmodel.GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := cmodel.F0()
	bs := cmodel.NewBranchSite(cali, t, cf, false)

	ds := optimize.NewMH(false, 0)
	ds.SetOptimizable(bs)
	ds.Run(5)
}

func TestAnnealing(tst *testing.T) {
	t, cali, err := cmodel.GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := cmodel.F3X4(cali)
	m0 := cmodel.NewM0(cali, t, cf)
	npar := len(m0.GetFloatParameters())
	if npar != 2 {
		tst.Error("Wrong number of parameters for M0:", npar)
	}

	m0.SetOptimizeBranchLengths()
	npar = len(m0.GetFloatParameters())
	if npar != 2+t.NNodes()-1 {
		tst.Error("Wrong number of parameters for M0 + branch lengths:", npar)
	}

	m0.SetParameters(2, 0.5)

	ds := optimize.NewMH(true, 0)
	ds.SetOptimizable(m0)
	ds.Run(5)
}

func TestLBFGSB(tst *testing.T) {
	t, cali, err := cmodel.GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := cmodel.F0()
	m0 := cmodel.NewM0(cali, t, cf)
	m0.SetParameters(2, 0.5)

	ds := optimize.NewLBFGSB()
	ds.SetOptimizable(m0)
	ds.Run(5)
}

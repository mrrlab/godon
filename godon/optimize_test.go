package main

import (
	"testing"

	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/codon"
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

	cf := codon.F3X4(cali)
	m0 := cmodel.NewM0(cali, t, cf)
	m0.SetParameters(2, 0.5)
	m := optimize.Optimizable(m0).Copy()
	npar := len(m0.GetFloatParameters())
	if npar != 2 {
		tst.Error("Wrong number of parameters for M0:", npar)
	}

	ds := optimize.NewDS()
	ds.SetOptimizable(m)
	ds.Quiet = true
	ds.Run(5)

	m = m.Copy()
	npar1 := len(m0.GetFloatParameters())
	npar2 := len(m.GetFloatParameters())
	if npar1 != npar2 {
		tst.Error("Parameter number mismatch after copy:", npar1, npar2)
	}

	ds = optimize.NewDS()
	ds.SetOptimizable(m)
	ds.Quiet = true
	ds.Run(5)
}

func TestMH(tst *testing.T) {
	t, cali, err := cmodel.GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := codon.F0(cali)
	bs := cmodel.NewBranchSite(cali, t, cf, false)

	mh := optimize.NewMH(false, 0)
	mh.SetOptimizable(bs)
	mh.Quiet = true
	mh.Run(5)

	m := bs.Copy()
	npar1 := len(bs.GetFloatParameters())
	npar2 := len(m.GetFloatParameters())
	if npar1 != npar2 {
		tst.Error("Parameter number mismatch after copy:", npar1, npar2)
	}

	mh = optimize.NewMH(false, 0)
	mh.SetOptimizable(m)
	mh.Quiet = true
	mh.Run(5)
}

func TestAnnealing(tst *testing.T) {
	t, cali, err := cmodel.GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := codon.F3X4(cali)
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

	an := optimize.NewMH(true, 0)
	an.SetOptimizable(m0)
	an.Quiet = true
	an.Run(5)

	m := m0.Copy()
	npar1 := len(m0.GetFloatParameters())
	npar2 := len(m.GetFloatParameters())
	if npar1 != npar2 {
		tst.Error("Parameter number mismatch after copy:", npar1, npar2)
	}
	an = optimize.NewMH(false, 0)
	an.SetOptimizable(m)
	an.Quiet = true
	an.Run(5)

}

func TestLBFGSB(tst *testing.T) {
	t, cali, err := cmodel.GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := codon.F0(cali)
	m0 := cmodel.NewM0(cali, t, cf)
	m0.SetParameters(2, 0.5)

	l := optimize.NewLBFGSB()
	l.SetOptimizable(m0)
	l.Quiet = true
	l.Run(5)

	m := m0.Copy()
	npar1 := len(m0.GetFloatParameters())
	npar2 := len(m.GetFloatParameters())
	if npar1 != npar2 {
		tst.Error("Parameter number mismatch after copy:", npar1, npar2)
	}

	l = optimize.NewLBFGSB()
	l.SetOptimizable(m)
	l.Quiet = true
	l.Run(5)
}

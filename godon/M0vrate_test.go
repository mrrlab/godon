package main

import (
	"testing"

	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
)

func TestM0vrate(tst *testing.T) {
	t, cali, err := cmodel.GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := codon.F3X4(cali)
	m0v := cmodel.NewM0vrate(cali, t, cf)
	m := optimize.Optimizable(m0v).Copy()
	npar := len(m0v.GetFloatParameters())
	if npar != 3 {
		tst.Error("Wrong number of parameters for M0vrate:", npar)
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

	ds = optimize.NewDS()
	ds.SetOptimizable(m)
	ds.Quiet = true
	ds.Run(5)
}

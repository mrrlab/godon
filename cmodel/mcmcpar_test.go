package cmodel

import (
	"math"
	"testing"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
)

func TestBranchSiteReprM0D1(tst *testing.T) {
	t, cali, err := GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := codon.F3X4(cali)

	m0 := NewM0(cali, t, cf)
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(m0)
	chain.Quiet = true
	chain.Run(10)
	L := m0.Likelihood()

	// Now reproduce
	k, w := m0.GetParameters()

	t, cali, err = GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}
	cf = codon.F3X4(cali)
	m0 = NewM0(cali, t, cf)
	m0.SetParameters(k, w)
	newL := m0.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

func TestBranchSiteReprM0D2(tst *testing.T) {
	t, cali, err := GetTreeAlignment(data2)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := codon.F0(cali)

	m0 := NewM0(cali, t, cf)
	as := optimize.NewAdaptiveSettings()
	m0.SetAdaptive(as)
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(m0)
	chain.Quiet = true
	chain.Run(10)
	L := m0.Likelihood()

	// Now reproduce
	k, w := m0.GetParameters()

	t, cali, err = GetTreeAlignment(data2)
	if err != nil {
		tst.Error("Error: ", err)
	}
	cf = codon.F0(cali)
	m0 = NewM0(cali, t, cf)
	m0.SetParameters(k, w)
	newL := m0.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

func TestBranchSiteReprM0D3(tst *testing.T) {
	t, cali, err := GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := codon.F3X4(cali)

	m0 := NewM0(cali, t, cf)
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(m0)
	chain.Quiet = true
	chain.Run(20)
	L := m0.Likelihood()

	// Now reproduce
	k, w := m0.GetParameters()

	_, cali, err = GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}
	cf = codon.F3X4(cali)
	m0 = NewM0(cali, t, cf)
	m0.SetParameters(k, w)
	newL := m0.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

func TestBranchSiteReprBSD1(tst *testing.T) {
	if testing.Short() {
		tst.Skip("skipping test in short mode.")
	}

	t, cali, err := GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := codon.F3X4(cali)

	h1 := NewBranchSite(cali, t, cf, false)
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(h1)
	chain.Quiet = true
	chain.Run(10)
	L := h1.Likelihood()

	// Now reproduce
	k, w0, w2, p0, p1 := h1.GetParameters()

	t, cali, err = GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}
	cf = codon.F3X4(cali)
	h1 = NewBranchSite(cali, t, cf, false)
	h1.SetParameters(k, w0, w2, p0, p1)
	newL := h1.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

func TestBranchSiteReprBSD2(tst *testing.T) {
	if testing.Short() {
		tst.Skip("skipping test in short mode.")
	}

	t, cali, err := GetTreeAlignment(data2)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := codon.F0(cali)

	h1 := NewBranchSite(cali, t, cf, false)
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(h1)
	chain.Quiet = true
	chain.Run(10)
	L := h1.Likelihood()

	// Now reproduce
	k, w0, w2, p0, p1 := h1.GetParameters()
	tst.Log("par=", k, w0, w2, p0, p1)

	t, cali, err = GetTreeAlignment(data2)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf = codon.F0(cali)
	h1 = NewBranchSite(cali, t, cf, false)
	h1.SetParameters(k, w0, w2, p0, p1)
	newL := h1.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}
func TestBranchSiteReprBSD3(tst *testing.T) {
	if testing.Short() {
		tst.Skip("skipping test in short mode.")
	}

	t, cali, err := GetTreeAlignment(data3)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := codon.F3X4(cali)

	h1 := NewBranchSite(cali, t, cf, true)
	as := optimize.NewAdaptiveSettings()
	h1.SetAdaptive(as)
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(h1)
	chain.Quiet = true
	chain.Run(10)
	L := h1.Likelihood()

	// Now reproduce
	k, w0, w2, p0, p1 := h1.GetParameters()
	tst.Log("par=", k, w0, w2, p0, p1)

	_, cali, err = GetTreeAlignment(data3)
	if err != nil {
		tst.Error("Error: ", err)
	}
	cf = codon.F3X4(cali)
	h1 = NewBranchSite(cali, t, cf, false)
	h1.SetParameters(k, w0, w2, p0, p1)

	newL := h1.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

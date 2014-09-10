package cmodel

import (
	"math"
	"testing"

	"bitbucket.com/Davydov/golh/optimize"
)

func TestBranchSiteReprM0D1(tst *testing.T) {
	t, cali, err := getTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F3X4(cali)

	m0 := NewM0(cali, t, cf, false)
	chain := optimize.NewMH(m0)
	chain.Quiet = true
	chain.Run(10)
	L := chain.L

	// Now reproduce
	k, w := m0.GetParameters()

	t, cali, err = getTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}
	cf = F3X4(cali)
	m0 = NewM0(cali, t, cf, false)
	m0.SetParameters(k, w)
	newL := m0.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

func TestBranchSiteReprM0D2(tst *testing.T) {
	t, cali, err := getTreeAlignment(data2)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F0()

	m0 := NewM0(cali, t, cf, false)
	as := optimize.NewAdaptiveSettings()
	m0.SetAdaptive(as)
	chain := optimize.NewMH(m0)
	chain.Quiet = true
	chain.Run(10)
	L := chain.L

	// Now reproduce
	k, w := m0.GetParameters()

	t, cali, err = getTreeAlignment(data2)
	if err != nil {
		tst.Error("Error: ", err)
	}
	cf = F0()
	m0 = NewM0(cali, t, cf, false)
	m0.SetParameters(k, w)
	newL := m0.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

func TestBranchSiteReprM0D3(tst *testing.T) {
	t, cali, err := getTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F3X4(cali)

	m0 := NewM0(cali, t, cf, true)
	chain := optimize.NewMH(m0)
	chain.Quiet = true
	chain.Run(20)
	L := chain.L

	// Now reproduce
	k, w := m0.GetParameters()

	_, cali, err = getTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}
	cf = F3X4(cali)
	m0 = NewM0(cali, t, cf, false)
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

	t, cali, err := getTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F3X4(cali)

	h1 := NewBranchSite(cali, t, cf, false)
	chain := optimize.NewMH(h1)
	chain.Quiet = true
	chain.Run(10)
	L := chain.L

	// Now reproduce
	k, w0, w2, p0, p1 := h1.GetParameters()

	t, cali, err = getTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}
	cf = F3X4(cali)
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

	t, cali, err := getTreeAlignment(data2)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F0()

	h1 := NewBranchSite(cali, t, cf, false)
	chain := optimize.NewMH(h1)
	chain.Quiet = true
	chain.Run(10)
	L := chain.L

	// Now reproduce
	k, w0, w2, p0, p1 := h1.GetParameters()
	tst.Log("par=", k, w0, w2, p0, p1)

	t, cali, err = getTreeAlignment(data2)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf = F0()
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

	t, cali, err := getTreeAlignment(data3)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F3X4(cali)

	h1 := NewBranchSite(cali, t, cf, true)
	as := optimize.NewAdaptiveSettings()
	h1.SetAdaptive(as)
	chain := optimize.NewMH(h1)
	chain.Quiet = true
	chain.Run(10)
	L := chain.L

	// Now reproduce
	k, w0, w2, p0, p1 := h1.GetParameters()
	tst.Log("par=", k, w0, w2, p0, p1)

	_, cali, err = getTreeAlignment(data3)
	if err != nil {
		tst.Error("Error: ", err)
	}
	cf = F3X4(cali)
	h1 = NewBranchSite(cali, t, cf, false)
	h1.SetParameters(k, w0, w2, p0, p1)

	newL := h1.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

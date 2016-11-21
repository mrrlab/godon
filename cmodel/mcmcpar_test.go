package cmodel

import (
	"math"
	"testing"

	"bitbucket.org/Davydov/godon/optimize"
)

func TestBranchSiteReprM0D1(tst *testing.T) {
	data, err := GetTreeAlignment(data1, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	m0 := NewM0(data)
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(m0)
	chain.Quiet = true
	chain.Run(10)
	L := m0.Likelihood()

	// Now reproduce
	k, w := m0.GetParameters()

	data, err = GetTreeAlignment(data1, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	m0 = NewM0(data)
	m0.SetParameters(k, w)
	newL := m0.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

func TestBranchSiteReprM0D2(tst *testing.T) {
	data, err := GetTreeAlignment(data2, "F0")
	if err != nil {
		tst.Error("Error: ", err)
	}

	m0 := NewM0(data)
	as := optimize.NewAdaptiveSettings()
	m0.SetAdaptive(as)
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(m0)
	chain.Quiet = true
	chain.Run(10)
	L := m0.Likelihood()

	// Now reproduce
	k, w := m0.GetParameters()

	data, err = GetTreeAlignment(data2, "F0")
	if err != nil {
		tst.Error("Error: ", err)
	}
	m0 = NewM0(data)
	m0.SetParameters(k, w)
	newL := m0.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

func TestBranchSiteReprM0D3(tst *testing.T) {
	data, err := GetTreeAlignment(data1, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	m0 := NewM0(data)
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(m0)
	chain.Quiet = true
	chain.Run(20)
	L := m0.Likelihood()

	// Now reproduce
	k, w := m0.GetParameters()

	dataNew, err := GetTreeAlignment(data1, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}
	dataNew.Tree = data.Tree
	m0 = NewM0(dataNew)
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

	data, err := GetTreeAlignment(data1, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	h1 := NewBranchSite(data, false)
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(h1)
	chain.Quiet = true
	chain.Run(10)
	L := h1.Likelihood()

	// Now reproduce
	k, w0, w2, p0, p1 := h1.GetParameters()

	data, err = GetTreeAlignment(data1, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	h1 = NewBranchSite(data, false)
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

	data, err := GetTreeAlignment(data2, "F0")
	if err != nil {
		tst.Error("Error: ", err)
	}

	h1 := NewBranchSite(data, false)
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(h1)
	chain.Quiet = true
	chain.Run(10)
	L := h1.Likelihood()

	// Now reproduce
	k, w0, w2, p0, p1 := h1.GetParameters()
	tst.Log("par=", k, w0, w2, p0, p1)

	data, err = GetTreeAlignment(data2, "F0")
	if err != nil {
		tst.Error("Error: ", err)
	}

	h1 = NewBranchSite(data, false)
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

	data, err := GetTreeAlignment(data3, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	h1 := NewBranchSite(data, true)
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

	dataNew, err := GetTreeAlignment(data3, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}
	dataNew.Tree = data.Tree

	h1 = NewBranchSite(dataNew, false)
	h1.SetParameters(k, w0, w2, p0, p1)

	newL := h1.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

package main

import (
	"math"
	"testing"

	"bitbucket.com/Davydov/golh/mcmc"
)

func TestBranchSiteReprBSD1(tst *testing.T) {
	if testing.Short() {
		tst.Skip("skipping test in short mode.")
	}

	t, cali, err := GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F3X4(cali)

	h1 := NewBranchSite(cali, t, cf, false)
	chain := mcmc.NewMH(h1)
	chain.Quiet = true
	chain.Run(10)
	L := chain.L

	// Now reproduce
	k, w0, w2, p0, p1 := h1.GetParameters()

	t, cali, err = GetTreeAlignment(data1)
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

	t, cali, err := GetTreeAlignment(data2)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F0()

	h1 := NewBranchSite(cali, t, cf, false)
	chain := mcmc.NewMH(h1)
	chain.Quiet = true
	chain.Run(10)
	L := chain.L

	// Now reproduce
	k, w0, w2, p0, p1 := h1.GetParameters()
	tst.Log("par=", k, w0, w2, p0, p1)

	t, cali, err = GetTreeAlignment(data2)
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

	t, cali, err := GetTreeAlignment(data3)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F3X4(cali)

	h1 := NewBranchSite(cali, t, cf, true)
	chain := mcmc.NewMH(h1)
	chain.Quiet = true
	chain.Run(10)
	L := chain.L

	// Now reproduce
	k, w0, w2, p0, p1 := h1.GetParameters()
	tst.Log("par=", k, w0, w2, p0, p1)

	_, cali, err = GetTreeAlignment(data3)
	if err != nil {
		tst.Error("Error: ", err)
	}
	cf = F3X4(cali)
	h1 = NewBranchSite(cali, t, cf, false)

	newL := h1.Likelihood()

	tst.Log("L=", newL, ", Ref=", L, ", diff=", math.Abs(L-newL))
	if math.IsNaN(L) || math.Abs(L-newL) > smallDiff {
		tst.Error("Expected ", L, ", got", newL)
	}
}

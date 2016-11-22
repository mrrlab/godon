package cmodel

import (
	"math"
	"testing"

	"github.com/op/go-logging"
)

const (
	data1 = "EMGT00050000008747.Drosophila.002"
	data2 = "ENSGT00550000073950.Euteleostomi.07.001"
	data3 = "EMGT00050000000025.Drosophila.001"

	// smallDiff is a threshold for testing
	// if likelihood ratio is larger error is emmited
	smallDiff = 1e-3
)

func init() {
	// disable logging for benchmarks
	logging.SetLevel(logging.WARNING, "optimize")
	logging.SetLevel(logging.ERROR, "cmodel")
}

/*** Test M0 ***/
func TestM0F0D1(tst *testing.T) {
	data, err := GetTreeAlignment(data1, "F0")
	if err != nil {
		tst.Error("Error: ", err)
	}

	m0 := NewM0(data)
	m0.SetParameters(2, 0.5)

	L := m0.Likelihood()
	refL := -2836.196647

	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.IsNaN(L) || math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestM0F3X4D1(tst *testing.T) {
	data, err := GetTreeAlignment(data1, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	m0 := NewM0(data)
	m0.SetParameters(2, 0.5)

	L := m0.Likelihood()
	refL := -2892.446106

	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.IsNaN(L) || math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestM0F0D2(tst *testing.T) {
	data, err := GetTreeAlignment(data2, "F0")
	if err != nil {
		tst.Error("Error: ", err)
	}

	m0 := NewM0(data)
	m0.SetParameters(1.79668, 0.09879)

	L := m0.Likelihood()
	refL := -1463.253413

	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.IsNaN(L) || math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestM0F3X4D2(tst *testing.T) {
	data, err := GetTreeAlignment(data2, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	m0 := NewM0(data)
	m0.SetParameters(3.12566, 0.03430)

	L := m0.Likelihood()
	refL := -1473.371833

	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.IsNaN(L) || math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestM0F0D3(tst *testing.T) {
	if testing.Short() {
		tst.Skip("skipping test in short mode.")
	}
	data, err := GetTreeAlignment(data3, "F0")
	if err != nil {
		tst.Error("Error: ", err)
	}

	m0 := NewM0(data)
	m0.SetParameters(1.77621, 0.10313)

	L := m0.Likelihood()
	refL := -48631.160712

	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.IsNaN(L) || math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

/*** Test BranchSite ***/
func TestBranchSiteF0D1(tst *testing.T) {
	data, err := GetTreeAlignment(data1, "F0")
	if err != nil {
		tst.Error("Error: ", err)
	}

	p0, p1 := 0.946800, 0.000098

	h1 := NewBranchSite(data, false)
	h1.SetParameters(1.909912, 0.020004, 1.000000, p0, p1)
	L := h1.Likelihood()
	refL := -2467.931313

	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.IsNaN(L) || math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestBranchSiteF3X4D1(tst *testing.T) {
	data, err := GetTreeAlignment(data1, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	p0, p1 := 0.934681, 0.000000

	h1 := NewBranchSite(data, false)
	h1.SetParameters(1.572572, 0.015173, 1.000000, p0, p1)
	L := h1.Likelihood()

	refL := -2474.003708
	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.IsNaN(L) || math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestBranchSiteF0D2(tst *testing.T) {
	data, err := GetTreeAlignment(data2, "F0")
	if err != nil {
		tst.Error("Error: ", err)
	}

	p0, p1 := 0.899776, 0.041502

	h1 := NewBranchSite(data, false)
	h1.SetParameters(1.876889, 0.014541, 1.682486, p0, p1)
	L := h1.Likelihood()
	refL := -1405.850679

	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.IsNaN(L) || math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestBranchSiteF3X4D2(tst *testing.T) {
	data, err := GetTreeAlignment(data2, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	p0, p1 := 0.932349, 0.005383

	h1 := NewBranchSite(data, false)
	h1.SetParameters(3.462166, 0.017381, 1.000000, p0, p1)
	refL := -1428.862660
	L := h1.Likelihood()

	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.IsNaN(L) || math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestBranchSiteF0D3(tst *testing.T) {
	if testing.Short() {
		tst.Skip("skipping test in short mode.")
	}
	data, err := GetTreeAlignment(data3, "F0")
	if err != nil {
		tst.Error("Error: ", err)
	}

	p0, p1 := 0.883725, 0.099870

	h1 := NewBranchSite(data, false)
	h1.SetParameters(1.94922, 0.072076, 13.936497, p0, p1)
	L := h1.Likelihood()
	refL := -48019.677814

	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.IsNaN(L) || math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestBranchSiteF3X4D3(tst *testing.T) {
	if testing.Short() {
		tst.Skip("skipping test in short mode.")
	}
	data, err := GetTreeAlignment(data3, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}

	p0, p1 := 0.902686, 0.079667

	h1 := NewBranchSite(data, false)
	h1.SetParameters(1.868810, 0.064487, 6.972271, p0, p1)
	refL := -47904.801502
	L := h1.Likelihood()

	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.IsNaN(L) || math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

/*** Benchmark M0 ***/
func BenchmarkM0F0D1(b *testing.B) {
	data, err := GetTreeAlignment(data1, "F0")
	if err != nil {
		b.Error("Error: ", err)
	}

	m0 := NewM0(data)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		m0.SetParameters(2, 0.5)
		m0.Likelihood()
	}
}

func BenchmarkM0F0D2(b *testing.B) {
	data, err := GetTreeAlignment(data2, "F0")
	if err != nil {
		b.Error("Error: ", err)
	}

	m0 := NewM0(data)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		m0.SetParameters(1.77621, 0.10313)
		m0.Likelihood()
	}
}

/*** Benchmark BranchSite ***/
func BenchmarkBranchSiteF0D1(b *testing.B) {
	data, err := GetTreeAlignment(data1, "F0")
	if err != nil {
		b.Error("Error: ", err)
	}

	p0, p1 := 0.946800, 0.000098

	h1 := NewBranchSite(data, false)
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		h1.SetParameters(1.909912, 0.020004, 1.000000, p0, p1)
		h1.Likelihood()
	}
}

func BenchmarkBranchSiteF0D2(b *testing.B) {
	data, err := GetTreeAlignment(data1, "F0")
	if err != nil {
		b.Error("Error: ", err)
	}

	p0, p1 := 0.899776, 0.041502

	h1 := NewBranchSite(data, false)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		h1.SetParameters(1.876889, 0.014541, 1.682486, p0, p1)
		h1.Likelihood()
	}
}

package main

import (
	"os"
	"math"
	"testing"
	"path"

	"bitbucket.com/Davydov/golh/bio"
	"bitbucket.com/Davydov/golh/tree"
)

const (
	data1 = "EMGT00050000008747.Drosophila.002"
	data2 = "ENSGT00550000073950.Euteleostomi.07.001"
	data3 = "EMGT00050000000025.Drosophila.001"

	smallDiff = 1e-4
)

func init() {
}

func GetTreeAlignment(data string) (t *tree.Tree, cali CodonSequences, err error) {
	tf, err := os.Open(path.Join("testdata", data +".nwk"))
	if err != nil {
		return
	}
	defer tf.Close()

	t, err = tree.ParseNewick(tf)
	if err != nil {
		return
	}

	af, err := os.Open(path.Join("testdata", data +".fst"))
	if err != nil {
		return
	}
	defer af.Close()

	ali, err := bio.ParseFasta(af)
	if err != nil {
		return
	}

	nCodon = initCodon()

	cali, err = ToCodonSequences(ali)
	if err != nil {
		return
	}
	nm2id = make(map[string]int)
	for i, s := range cali {
		nm2id[s.Name] = i
	}

	return
}

func TestM0_1(tst *testing.T) {
	t, cali, err := GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F0()
	L := M0(cali, t, cf, 2, 0.5)
	refL := -2836.196647
	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestM0_2(tst *testing.T) {
	t, cali, err := GetTreeAlignment(data2)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F0()
	L := M0(cali, t, cf, 1.79668, 0.09879)
	refL := -1463.253413
	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestM0_3(tst *testing.T) {
	if testing.Short() {
		tst.Skip("skipping test in short mode.")
	}
	t, cali, err := GetTreeAlignment(data3)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F0()
	L := M0(cali, t, cf, 1.77621, 0.10313)
	refL := -48631.160712
	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestH1_1(tst *testing.T) {
	t, cali, err := GetTreeAlignment(data1)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F0()
	p0, p1 := 0.946800, 0.000098
	p2a := (1 - p0 - p1) * p0 / (p0 + p1)
	p2b := (1 - p0 - p1) * p1 / (p0 + p1)
	L := H1(cali, t, cf, 1.909912, 0.020004, 1.000000, p0, p1, p2a, p2b)
	refL := -2467.931313
	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestH1_2(tst *testing.T) {
	t, cali, err := GetTreeAlignment(data2)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F0()
	p0, p1 := 0.899776, 0.041502
	p2a := (1 - p0 - p1) * p0 / (p0 + p1)
	p2b := (1 - p0 - p1) * p1 / (p0 + p1)
	L := H1(cali, t, cf, 1.876889, 0.014541, 1.682486, p0, p1, p2a, p2b)
	refL := -1405.850679
	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func TestH1_3(tst *testing.T) {
	if testing.Short() {
		tst.Skip("skipping test in short mode.")
	}
	t, cali, err := GetTreeAlignment(data3)
	if err != nil {
		tst.Error("Error: ", err)
	}

	cf := F0()
	p0, p1 := 0.883725, 0.099870
	p2a := (1 - p0 - p1) * p0 / (p0 + p1)
	p2b := (1 - p0 - p1) * p1 / (p0 + p1)
	L := H1(cali, t, cf, 1.94922, 0.072076, 13.936497, p0, p1, p2a, p2b)
	refL := -48019.677814
	tst.Log("L=", L, ", Ref=", refL, ", diff=", math.Abs(L-refL))
	if math.Abs(L-refL) > smallDiff {
		tst.Error("Expected ", refL, ", got", L)
	}
}

func BenchmarkM0_1(b *testing.B) {
	t, cali, err := GetTreeAlignment(data1)
	if err != nil {
		b.Error("Error: ", err)
	}

	cf := F0()

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		M0(cali, t, cf, 2, 0.5)
	}
}

func BenchmarkM0_2(b *testing.B) {
	t, cali, err := GetTreeAlignment(data2)
	if err != nil {
		b.Error("Error: ", err)
	}

	cf := F0()

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		M0(cali, t, cf, 1.77621, 0.10313)
	}
}

func BenchmarkH1_1(b *testing.B) {
	t, cali, err := GetTreeAlignment(data1)
	if err != nil {
		b.Error("Error: ", err)
	}

	cf := F0()
	p0, p1 := 0.946800, 0.000098
	p2a := (1 - p0 - p1) * p0 / (p0 + p1)
	p2b := (1 - p0 - p1) * p1 / (p0 + p1)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		H1(cali, t, cf, 1.909912, 0.020004, 1.000000, p0, p1, p2a, p2b)
	}
}

func BenchmarkH1_2(b *testing.B) {
	t, cali, err := GetTreeAlignment(data1)
	if err != nil {
		b.Error("Error: ", err)
	}

	cf := F0()
	p0, p1 := 0.899776, 0.041502
	p2a := (1 - p0 - p1) * p0 / (p0 + p1)
	p2b := (1 - p0 - p1) * p1 / (p0 + p1)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		H1(cali, t, cf, 1.876889, 0.014541, 1.682486, p0, p1, p2a, p2b)
	}
}

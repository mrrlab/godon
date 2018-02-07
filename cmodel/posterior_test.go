package cmodel

import (
	"math"
	"testing"
)

func comparePosterior(posterior []float64, reference map[int]float64, tst *testing.T) {
	for i, p := range posterior {
		r, ok := reference[i+1]
		if ok {
			diff := math.Abs(p - r)
			tst.Log("i=", i, "p=", p, "r=", r, "diff=", diff)
			if diff > smallDiff {
				tst.Error("Expected ", r, ", got", p)
			}
		} else {
			diff := r - 0.5
			if diff > smallDiff {
				tst.Error("Expected <0.5, got", p)
			}
		}
	}
}

func TestBranchSiteNebD3(tst *testing.T) {
	data, err := GetTreeAlignment(data3, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}
	h1 := NewBranchSite(data, false)
	h1.SetParameters(1.86881, 0.06449, 6.97229, 0.90269, 0.07967)
	h1.Final(true, false, false, false, false)
	neb := h1.summary.SitePosteriorNEB
	referenceNeb := map[int]float64{
		// this comes from PAML
		300:  0.934,
		315:  0.950,
		361:  0.555,
		487:  0.527,
		670:  0.964,
		711:  0.635,
		718:  0.962,
		735:  0.842,
		736:  0.787,
		816:  0.815,
		1017: 0.686,
		1052: 0.886,
		1079: 0.760,
		1100: 0.777,
		1302: 0.982,
		1315: 0.925,
		1399: 0.952,
		1656: 0.902,
		1906: 0.635,
		2138: 0.662,
		2322: 0.785,
		2400: 0.761,
		2419: 0.968,
		2747: 0.878,
		2758: 0.932,
	}
	comparePosterior(neb, referenceNeb, tst)

}

func TestBranchSiteBebD3(tst *testing.T) {
	if testing.Short() {
		tst.Skip("skipping test in short mode.")
	}

	data, err := GetTreeAlignment(data3, "F3X4")
	if err != nil {
		tst.Error("Error: ", err)
	}
	h1 := NewBranchSite(data, false)
	h1.SetParameters(1.86881, 0.06449, 6.97229, 0.90269, 0.07967)
	h1.Final(false, true, false, false, false)
	beb := h1.summary.SitePosteriorBEB
	referenceBeb := map[int]float64{
		// this comes from PAML
		217:  0.557,
		300:  0.971,
		315:  0.971,
		357:  0.550,
		361:  0.720,
		670:  0.972,
		711:  0.785,
		718:  0.969,
		723:  0.518,
		735:  0.858,
		736:  0.835,
		816:  0.637,
		1017: 0.788,
		1051: 0.506,
		1052: 0.763,
		1079: 0.630,
		1100: 0.713,
		1302: 0.988,
		1315: 0.934,
		1399: 0.977,
		1534: 0.643,
		1656: 0.852,
		1711: 0.503,
		1906: 0.666,
		2131: 0.567,
		2138: 0.663,
		2322: 0.800,
		2400: 0.830,
		2418: 0.611,
		2419: 0.940,
		2461: 0.606,
		2747: 0.948,
		2758: 0.953,
	}
	comparePosterior(beb, referenceBeb, tst)
}

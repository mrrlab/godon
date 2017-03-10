package cmodel

// This is test for highly variable sequence alignment. data4 was
// generates using cosim:
//  cosim --model M0 --omega 0.2 --kappa 2 --codon-alpha .1 t.br.nwk 500 ali.fst
// This variability between positions gives zero probabiltiy for certain
// positions when branches are short. Zero probability produce -Inf lnL.

import (
	"math"
	"testing"
)

const data4 = "m0g"

/*** Test M0 ***/
func TestM0GSim(tst *testing.T) {
	data, err := GetTreeAlignment(data4, "F0")
	if err != nil {
		tst.Error("Error: ", err)
	}

	m0 := NewM0(data)
	m0.SetParameters(1, 1)

	L := m0.Likelihood()
	if math.IsNaN(L) || math.IsInf(L, 0) || L > 0 {
		tst.Error("Invalid likelihood:", L)
	}
}

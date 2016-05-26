package cmodel

import (
	"testing"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
)

func BenchmarkMCMCD1(b *testing.B) {
	t, cali, err := GetTreeAlignment(data1)
	if err != nil {
		b.Error("Error: ", err)
	}

	cf := codon.F0()

	m0 := NewM0(cali, t, cf)

	b.ResetTimer()

	m0.SetDefaults()
	chain := optimize.NewMH(false, 0)
	chain.SetOptimizable(m0)
	chain.Quiet = true
	chain.Run(100)
}

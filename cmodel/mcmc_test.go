package cmodel

import (
	"testing"

	"bitbucket.org/Davydov/godon/optimize"
)

func BenchmarkMCMCD1(b *testing.B) {
	data, err := GetTreeAlignment(data1, "F0")
	if err != nil {
		b.Error("Error: ", err)
	}

	m0 := NewM0(data)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		m0.SetDefaults()
		chain := optimize.NewMH(false, 0)
		chain.SetOptimizable(m0)
		chain.Quiet = true
		chain.Run(100)
	}
}

package main

import (
	"testing"

	"bitbucket.com/Davydov/golh/mcmc"
)

func BenchmarkMCMCD1(b *testing.B) {
	t, cali, err := GetTreeAlignment(data1)
	if err != nil {
		b.Error("Error: ", err)
	}

	cf := F0()

	m0 := NewM0(cali, t, cf)

	b.ResetTimer()

	m0.SetDefaults()
	mcmc.MCMC(m0, 100, 0, 1000)
}


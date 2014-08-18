package main

import (
	"log"
	"math"
	"math/rand"
	"time"
)

const stdev = 1e-2

func init() {
	rand.Seed(time.Now().UnixNano())
}

func MCMC(m Optimizable, iterations int) {
	accepted := 0
	np := m.GetNumberOfParameters()
	L := m.Likelihood()
	for i := 0; i < iterations; i++ {
		p := rand.Intn(np)
		val := m.GetParameter(p)
		newVal := val + rand.NormFloat64()*stdev
		m.SetParameter(p, newVal)
		newL := m.Likelihood()
		a := math.Exp(newL - L)
		if a < 1 && rand.Float64() > a {
			m.SetParameter(p, val)
		} else {
			L = newL
			accepted++
			log.Printf("%d: Accept, L=%f", i, newL)
		}
	}
	log.Printf("Finished MCMC, acceptance rate %f%%", 100*float64(accepted)/float64(iterations))
	for i := 0; i < np; i++ {
		log.Printf("%s=%f", m.GetParameterName(i), m.GetParameter(i))
	}
}

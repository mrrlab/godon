package main

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"strconv"
	"time"
)

const stdev = 1e-2

func init() {
	rand.Seed(time.Now().UnixNano())
}

func MCMC(m Optimizable, burnIn, iterations int, report int) {
	accepted := 0
	np := m.GetNumberOfParameters()
	L := m.Likelihood()
	if burnIn > 0 {
		log.Printf("Burnin for %d iterations", burnIn)
	}
	fmt.Printf("iteration likelihood %s\n", ParameterNamesString(m))
	for i := 0; i < burnIn+iterations+burnIn; i++ {
		iter := i - burnIn
		if iter == 0 {
			log.Print("Starting sampling")
		}
		if iter >= 0 && iter%report == 0 {
			log.Printf("%d: L=%f", i, L)
			fmt.Printf("%d %f %s\n", i, L, ParameterString(m))
		}
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
		}
	}
	log.Printf("Finished MCMC, acceptance rate %f%%", 100*float64(accepted)/float64(iterations+burnIn))
	for i := 0; i < np; i++ {
		log.Printf("%s=%f", m.GetParameterName(i), m.GetParameter(i))
	}
}

func ParameterNamesString(m Optimizable) (s string) {
	for i := 0; i < m.GetNumberOfParameters(); i++ {
		if i != 0 {
			s += " "
		}
		s += m.GetParameterName(i)
	}
	return
}
func ParameterString(m Optimizable) (s string) {
	for i := 0; i < m.GetNumberOfParameters(); i++ {
		if i != 0 {
			s += " "
		}
		s += strconv.FormatFloat(m.GetParameter(i), 'f', 6, 64)
	}
	return
}

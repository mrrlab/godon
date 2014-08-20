package mcmc

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"strconv"
	"time"
)

const (
	STDEV = 1e-2
)

type Optimizable interface {
	SetDefaults()
	GetNumberOfParameters() int
	GetParameterName(i int) string
	GetParameter(i int) float64
	SetParameter(i int, val float64)
	Likelihood() float64
}

func init() {
	rand.Seed(time.Now().UnixNano())
}

func MCMC(m Optimizable, burnIn, iterations int, report int) {
	np := m.GetNumberOfParameters()
	L := m.Likelihood()
	if burnIn > 0 {
		log.Printf("Burnin for %d iterations", burnIn)
	}
	accepted := 0
	fmt.Printf("iteration\tlikelihood\t%s\n",ParameterNamesString(m))
	for i := 0; i < burnIn+iterations+burnIn; i++ {
		if i%report == 0 && i > 0 {
			log.Printf("Acceptance rate %f%%", 100*float64(accepted)/float64(report))
			accepted = 0
		}

		iter := i - burnIn
		if iter == 0 {
			log.Print("Starting sampling")
		}
		if iter >= 0 && iter%report == 0 {
			log.Printf("%d: L=%f", i, L)
			fmt.Printf("%d\t%f\t%s\n", i, L, ParameterString(m))
		}
		p := rand.Intn(np)
		val := m.GetParameter(p)
		newVal := val + rand.NormFloat64()*STDEV
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
	log.Print("Finished MCMC")
	for i := 0; i < np; i++ {
		log.Printf("%s=%f", m.GetParameterName(i), m.GetParameter(i))
	}
}

func ParameterNamesString(m Optimizable) (s string) {
	for i := 0; i < m.GetNumberOfParameters(); i++ {
		if i != 0 {
			s += "\t"
		}
		s += m.GetParameterName(i)
	}
	return
}
func ParameterString(m Optimizable) (s string) {
	for i := 0; i < m.GetNumberOfParameters(); i++ {
		if i != 0 {
			s += "\t"
		}
		s += strconv.FormatFloat(m.GetParameter(i), 'f', 6, 64)
	}
	return
}

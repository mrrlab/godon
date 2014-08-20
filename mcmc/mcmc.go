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

type MCMC struct {
	Optimizable
	L float64
	np int
	i int
	RepPeriod int
	AccPeriod int
}

func NewMCMC(m Optimizable) *MCMC {
	return &MCMC{Optimizable: m,
	RepPeriod: 1000,
	np:  m.GetNumberOfParameters(),
	AccPeriod: 1000,
	}
}

func (m *MCMC) Run(iterations int) {
	m.L = m.Likelihood()
	m.PrintHeader()
	accepted := 0
	for m.i = 0; m.i < iterations; m.i++ {
		if m.i > 0 && m.i%m.AccPeriod == 0 {
			log.Printf("Acceptance rate %f%%", 100*float64(accepted)/float64(m.RepPeriod))
			accepted = 0
		}

		if m.i%m.RepPeriod == 0 {
			log.Printf("%d: L=%f", m.i, m.L)
			m.PrintLine()
		}
		p := rand.Intn(m.np)
		val := m.GetParameter(p)
		newVal := val + rand.NormFloat64()*STDEV
		m.SetParameter(p, newVal)
		newL := m.Likelihood()
		a := math.Exp(newL - m.L)
		if a < 1 && rand.Float64() > a {
			m.SetParameter(p, val)
		} else {
			m.L = newL
			accepted++
		}
	}
	log.Print("Finished MCMC")
	m.PrintFinal()
}

func (m *MCMC) PrintHeader() {
	fmt.Printf("iteration\tlikelihood\t%s\n",ParameterNamesString(m))
}

func (m *MCMC) PrintLine() {
	fmt.Printf("%d\t%f\t%s\n", m.i, m.L, ParameterString(m))
}

func (m *MCMC) PrintFinal() {
	for i := 0; i < m.np; i++ {
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

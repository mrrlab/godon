package mcmc

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"strconv"
	"time"
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
	L         float64
	np        int
	i         int
	RepPeriod int
	AccPeriod int
	SD        float64
	pnames    []string
	*Adaptive
}

func NewMCMC(m Optimizable) (mcmc *MCMC) {
	np := m.GetNumberOfParameters()
	mcmc = &MCMC{Optimizable: m,
		np:        np,
		RepPeriod: 10,
		AccPeriod: 10,
		SD:        1e-2,
		pnames:    make([]string, np),
	}
	for p := 0; p < np; p++ {
		mcmc.pnames[p] = m.GetParameterName(p)
	}
	return
}

func (m *MCMC) SetAdaptive(adaptive bool) {
	if adaptive {
		log.Print("Setting adaptive")
		m.Adaptive = NewAdaptive(m.np, m.pnames, m.SD)
	} else {
		log.Print("Setting nonadaptive")
		m.Adaptive = nil
	}
}

func (m *MCMC) Run(iterations int) {
	m.L = m.Likelihood()
	m.PrintHeader()
	accepted := 0
	for m.i = 0; m.i < iterations; m.i++ {
		if m.i > 0 && m.i%m.AccPeriod == 0 {
			log.Printf("Acceptance rate %.2f%%", 100*float64(accepted)/float64(m.AccPeriod))
			accepted = 0
		}

		if m.i%m.RepPeriod == 0 {
			log.Printf("%d: L=%f", m.i, m.L)
			m.PrintLine()
		}
		p := rand.Intn(m.np)
		val := m.GetParameter(p)
		var sd float64
		if m.Adaptive != nil {
			sd = m.Adaptive.SD[p]
		} else {
			sd = m.SD
		}
		newVal := val + rand.NormFloat64()*sd
		m.SetParameter(p, newVal)
		newL := m.Likelihood()
		a := math.Exp(newL - m.L)
		if a < 1 && rand.Float64() > a {
			m.SetParameter(p, val)
		} else {
			if m.Adaptive != nil {
				m.UpdateMu(p, newVal)
			}
			m.L = newL
			accepted++
		}
	}
	log.Print("Finished MCMC")
	m.PrintFinal()
}

func (m *MCMC) PrintHeader() {
	fmt.Printf("iteration\tlikelihood\t%s\n", m.ParameterNamesString())
}

func (m *MCMC) PrintLine() {
	fmt.Printf("%d\t%f\t%s\n", m.i, m.L, m.ParameterString())
}

func (m *MCMC) PrintFinal() {
	for i := 0; i < m.np; i++ {
		log.Printf("%s=%f", m.pnames[i], m.GetParameter(i))
	}
}

func (m *MCMC) ParameterNamesString() (s string) {
	for i := 0; i < m.np; i++ {
		if i != 0 {
			s += "\t"
		}
		s += m.pnames[i]
	}
	return
}

func (m *MCMC) ParameterString() (s string) {
	for i := 0; i < m.np; i++ {
		if i != 0 {
			s += "\t"
		}
		s += strconv.FormatFloat(m.GetParameter(i), 'f', 6, 64)
	}
	return
}

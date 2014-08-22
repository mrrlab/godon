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

type MH struct {
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

func NewMH(m Optimizable) (mcmc *MH) {
	np := m.GetNumberOfParameters()
	mcmc = &MH{Optimizable: m,
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

func (m *MH) SetAdaptive(ap *AdaptiveParameters) {
	if ap != nil {
		log.Print("Setting adaptive")
		m.Adaptive = NewAdaptive(m.np, m.pnames, m.SD, ap)
	} else {
		log.Print("Setting nonadaptive")
	}
}

func (m *MH) Run(iterations int) {
	m.L = m.Likelihood()
	if m.Adaptive != nil {
		log.Print(m.Adaptive)
	}
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
			if m.Adaptive != nil && m.i >= m.Skip {
				m.UpdateMu(p, newVal)
			}
			m.L = newL
			accepted++
		}
	}
	log.Print("Finished MCMC")
	m.PrintFinal()
}

func (m *MH) PrintHeader() {
	fmt.Printf("iteration\tlikelihood\t%s\n", m.ParameterNamesString())
}

func (m *MH) PrintLine() {
	fmt.Printf("%d\t%f\t%s\n", m.i, m.L, m.ParameterString())
}

func (m *MH) PrintFinal() {
	for i := 0; i < m.np; i++ {
		log.Printf("%s=%f", m.pnames[i], m.GetParameter(i))
	}
}

func (m *MH) ParameterNamesString() (s string) {
	for i := 0; i < m.np; i++ {
		if i != 0 {
			s += "\t"
		}
		s += m.pnames[i]
	}
	return
}

func (m *MH) ParameterString() (s string) {
	for i := 0; i < m.np; i++ {
		if i != 0 {
			s += "\t"
		}
		s += strconv.FormatFloat(m.GetParameter(i), 'f', 6, 64)
	}
	return
}

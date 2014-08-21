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
	*Adaptive
}

func NewMCMC(m Optimizable) *MCMC {
	return &MCMC{Optimizable: m,
		np:        m.GetNumberOfParameters(),
		RepPeriod: 10,
		AccPeriod: 10,
		SD:        1e-2,
	}
}

func (m *MCMC) SetAdaptive(adaptive bool) {
	if adaptive {
		log.Print("Setting adaptive")
		m.Adaptive = NewAdaptive(m.np, m.SD)
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
			log.Printf("Acceptance rate %f%%", 100*float64(accepted)/float64(m.RepPeriod))
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
	fmt.Printf("iteration\tlikelihood\t%s\n", ParameterNamesString(m))
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

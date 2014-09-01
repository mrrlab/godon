package mcmc

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"os/signal"
)

type Optimizable interface {
	GetModelParameters() Parameters
	Likelihood() float64
}

type MH struct {
	Optimizable
	L          float64
	MaxL       float64
	i          int
	RepPeriod  int
	AccPeriod  int
	SD         float64
	sig        chan os.Signal
	parameters Parameters
	Quiet      bool
}

func NewMH(m Optimizable) (mcmc *MH) {
	mcmc = &MH{Optimizable: m,
		RepPeriod:  10,
		AccPeriod:  10,
		SD:         1e-2,
		parameters: m.GetModelParameters(),
	}
	return
}

func (m *MH) WatchSignals(sigs ...os.Signal) {
	m.sig = make(chan os.Signal, 1)
	signal.Notify(m.sig, sigs...)
}

func (m *MH) Run(iterations int) {
	m.L = m.Likelihood()
	m.MaxL = m.L
	m.PrintHeader()
	accepted := 0
Iter:
	for m.i = 0; m.i < iterations; m.i++ {
		if !m.Quiet && m.i > 0 && m.i%m.AccPeriod == 0 {
			log.Printf("Acceptance rate %.2f%%", 100*float64(accepted)/float64(m.AccPeriod))
			accepted = 0
		}

		if !m.Quiet && m.i%m.RepPeriod == 0 {
			log.Printf("%d: L=%f", m.i, m.L)
			m.PrintLine()
		}
		p := rand.Intn(len(m.parameters))
		par := m.parameters[p]
		par.Propose()
		newL := m.Likelihood()

		a := math.Exp(par.Prior() - par.OldPrior() + newL - m.L)
		if a > 1 || rand.Float64() < a {
			m.L = newL
			par.Accept(m.i)
			accepted++
			m.MaxL = math.Max(m.L, m.MaxL)
		} else {
			par.Reject()
		}

		select {
		case s := <-m.sig:
			log.Printf("Received signal %v, exiting.", s)
			break Iter
		default:
		}
	}
	if !m.Quiet {
		log.Print("Finished MCMC")
		log.Printf("Maximum likelihood: %v", m.MaxL)
	}
	m.PrintFinal()
}

func (m *MH) PrintHeader() {
	if !m.Quiet {
		fmt.Printf("iteration\tlikelihood\t%s\n", m.ParameterNamesString())
	}
}

func (m *MH) PrintLine() {
	if !m.Quiet {
		fmt.Printf("%d\t%f\t%s\n", m.i, m.L, m.ParameterString())
	}
}

func (m *MH) PrintFinal() {
	if !m.Quiet {
		for _, par := range m.parameters {
			log.Printf("%s=%s", par.Name(), par.Value())
		}
	}
}

func (m *MH) ParameterNamesString() (s string) {
	for i, par := range m.parameters {
		if i != 0 {
			s += "\t"
		}
		s += par.Name()
	}
	return
}

func (m *MH) ParameterString() (s string) {
	for i, par := range m.parameters {
		if i != 0 {
			s += "\t"
		}
		s += par.Value()
	}
	return
}

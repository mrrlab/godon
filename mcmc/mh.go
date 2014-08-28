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
	GetParameters() Parameters
	Likelihood() float64
}

type MH struct {
	Optimizable
	L          float64
	i          int
	RepPeriod  int
	AccPeriod  int
	SD         float64
	sig        chan os.Signal
	parameters Parameters
}

func NewMH(m Optimizable) (mcmc *MH) {
	mcmc = &MH{Optimizable: m,
		RepPeriod:  10,
		AccPeriod:  10,
		SD:         1e-2,
		parameters: m.GetParameters(),
	}
	return
}

func (m *MH) WatchSignals(sigs ...os.Signal) {
	m.sig = make(chan os.Signal, 1)
	signal.Notify(m.sig, sigs...)
}

func (m *MH) SetAdaptive() {
	panic("not implemented")
	/*if ap != nil {
		log.Print("Setting adaptive")
		//m.Adaptive = NewAdaptive(m.np, m.pnames, m.SD, ap)
	} else {
		log.Print("Setting nonadaptive")
	}*/
}

func (m *MH) Run(iterations int) {
	m.L = m.Likelihood()
	m.PrintHeader()
	accepted := 0
Iter:
	for m.i = 0; m.i < iterations; m.i++ {
		if m.i > 0 && m.i%m.AccPeriod == 0 {
			log.Printf("Acceptance rate %.2f%%", 100*float64(accepted)/float64(m.AccPeriod))
			accepted = 0
		}

		if m.i%m.RepPeriod == 0 {
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
			par.Accept()
			accepted++
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
	for _, par := range m.parameters {
		log.Printf("%s=%s", par.Name(), par.Value())
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

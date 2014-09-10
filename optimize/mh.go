package optimize

import (
	"log"
	"math"
	"math/rand"
)

type Optimizable interface {
	GetModelParameters() Parameters
	Likelihood() float64
}

type MH struct {
	*BaseOptimizer
	AccPeriod int
	SD        float64
}

func NewMH() (mcmc *MH) {
	mcmc = &MH{
		BaseOptimizer: &BaseOptimizer{
			repPeriod: 10,
		},
		AccPeriod: 10,
		SD:        1e-2,
	}
	return
}

func (m *MH) Run(iterations int) {
	m.l = m.Likelihood()
	m.maxL = m.l
	m.maxLPar = m.ParameterString()
	m.PrintHeader()
	accepted := 0
Iter:
	for m.i = 0; m.i < iterations; m.i++ {
		if !m.Quiet && m.i > 0 && m.i%m.AccPeriod == 0 {
			log.Printf("Acceptance rate %.2f%%", 100*float64(accepted)/float64(m.AccPeriod))
			accepted = 0
		}

		if !m.Quiet && m.i%m.repPeriod == 0 {
			log.Printf("%d: L=%f", m.i, m.l)
			m.PrintLine()
		}
		p := rand.Intn(len(m.parameters))
		par := m.parameters[p]
		par.Propose()
		newL := m.Likelihood()

		a := math.Exp(par.Prior() - par.OldPrior() + newL - m.l)
		if a > 1 || rand.Float64() < a {
			m.l = newL
			par.Accept(m.i)
			accepted++
			if m.l > m.maxL {
				m.maxL = m.l
				m.maxLPar = m.ParameterString()
			}
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
		log.Printf("Maximum likelihood: %v", m.maxL)
		log.Printf("Parameter  names: %v", m.ParameterNamesString())
		log.Printf("Parameter values: %v", m.maxLPar)
	}
	m.PrintFinal()
}

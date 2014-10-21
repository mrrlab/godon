package optimize

import (
	"log"
	"math"
	"math/rand"
)

type MH struct {
	BaseOptimizer
	parameters FloatParameters
	Optimizable
	AccPeriod int
	annealing bool
	SD        float64
}

func NewMH(annealing bool) (mcmc *MH) {
	mcmc = &MH{
		BaseOptimizer: BaseOptimizer{
			repPeriod: 10,
		},
		AccPeriod: 10,
		SD:        1e-2,
		annealing: annealing,
	}
	return
}

func (m *MH) SetOptimizable(opt Optimizable) {
	m.Optimizable = opt
	m.parameters = opt.GetFloatParameters()
}

func (m *MH) Run(iterations int) {
	m.l = m.Likelihood()
	m.maxL = m.l
	m.maxLPar = m.parameters.Values(m.maxLPar)
	m.PrintHeader(m.parameters)
	accepted := 0
Iter:
	for m.i = 0; m.i < iterations; m.i++ {
		var T float64
		if m.annealing {
			T = math.Pow(0.9, float64(m.i)/float64(iterations)*100)
		} else {
			T = 1
		}
		if !m.Quiet && m.i > 0 && m.i%m.AccPeriod == 0 {
			log.Printf("Acceptance rate %.2f%%", 100*float64(accepted)/float64(m.AccPeriod))
			accepted = 0
		}

		if !m.Quiet && m.i%m.repPeriod == 0 {
			if m.annealing {
				log.Printf("%d: L=%f, T=%f", m.i, m.l, T)
			} else {
				log.Printf("%d: L=%f", m.i, m.l)
			}
			m.PrintLine(m.parameters, m.l)
		}
		p := rand.Intn(len(m.parameters))
		par := m.parameters[p]
		par.Propose()
		newL := m.Likelihood()

		var a float64
		if m.annealing {
			a = math.Exp((newL - m.l) / T)
		} else {
			a = math.Exp((par.Prior() - par.OldPrior() + newL - m.l))
		}

		if a > 1 || rand.Float64() < a {
			m.l = newL
			par.Accept(m.i)
			accepted++
			if m.l > m.maxL {
				m.maxL = m.l
				m.maxLPar = m.parameters.Values(m.maxLPar)
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
		log.Printf("Parameter  names: %v", m.parameters.NamesString())
		log.Printf("Parameter values: %v", m.GetMaxLParameters())
	}
	m.PrintFinal(m.parameters)
}

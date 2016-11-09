package optimize

import (
	"math"
	"math/rand"
)

// MH is a Metropolis-Hastings sampler.
type MH struct {
	BaseOptimizer
	AccPeriod int
	annealing bool
	// iteration to skip before annealing
	annealingSkip int
	SD            float64
}

// NewMH creates a new MH sampler.
func NewMH(annealing bool, annealingSkip int) (mcmc *MH) {
	mcmc = &MH{
		BaseOptimizer: BaseOptimizer{
			repPeriod: 10,
		},
		AccPeriod:     10,
		SD:            1e-2,
		annealing:     annealing,
		annealingSkip: annealingSkip,
	}
	return
}

// Run starts sampling.
func (m *MH) Run(iterations int) {
	m.SaveStart()
	m.PrintHeader()
	accepted := 0
	lastReported := -1
Iter:
	for m.i = 0; m.i < iterations; m.i++ {
		var T float64
		if m.annealing && m.i >= m.annealingSkip {
			T = math.Pow(0.9, float64(m.i-m.annealingSkip)/float64(iterations-m.annealingSkip)*100)
		} else {
			T = 1
		}
		if m.i > 0 && m.i%m.AccPeriod == 0 {
			log.Infof("Acceptance rate %.2f%%", 100*float64(accepted)/float64(m.AccPeriod))
			accepted = 0
		}

		if m.i%m.repPeriod == 0 {
			if m.annealing {
				log.Debugf("%d: L=%f, T=%f", m.i, m.l, T)
			} else {
				log.Debugf("%d: L=%f", m.i, m.l)
			}
			m.PrintLine(m.parameters, m.l)
			lastReported = m.i
		}
		p := rand.Intn(len(m.parameters))
		par := m.parameters[p]
		par.Propose()
		newL := m.Likelihood()
		m.calls++

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
			log.Warningf("Received signal %v, exiting.", s)
			break Iter
		default:
		}
	}

	if m.i != lastReported {
		m.PrintLine(m.parameters, m.l)
	}
}

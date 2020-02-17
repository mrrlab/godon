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
	l := m.startL
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

		m.PrintLine(m.parameters, l, m.repPeriod)
		if m.i%m.repPeriod == 0 {
			if m.annealing {
				log.Debugf("%d: L=%f, T=%f", m.i, l, T)
			} else {
				log.Debugf("%d: L=%f", m.i, l)
			}
			lastReported = m.i
		}
		p := rand.Intn(len(m.parameters))
		par := m.parameters[p]
		par.Propose()
		newL := m.Likelihood()
		m.calls++

		var a float64
		if m.annealing {
			a = math.Exp((newL - l) / T)
		} else {
			a = math.Exp((par.Prior() - par.OldPrior() + newL - l))
		}

		if a > 1 || rand.Float64() < a {
			l = newL
			par.Accept(m.i)
			accepted++
			if l > m.maxL {
				m.maxL = l
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
		m.PrintLine(m.parameters, l, 1)
	}

	m.SaveCheckpoint(true)
	m.saveDeltaT()
}

package mcmc

import (
	"log"
	"math"
)

type Adaptive struct {
	np         int
	pnames     []string
	t          []int
	loct       []int
	bmu        []float64
	mu         []float64
	vals       []chan float64
	sum        []float64
	sumsq      []float64
	SD         []float64
	delta      []bool
	converged  []bool
	nconverged int

	// constants
	wSize      int
	K          int
	Skip       int
	MaxUpdate  int
	Epsilon    float64
	C          float64
	Nu         float64
	Lambda     float64
}

func NewAdaptive(np int, pnames []string, sd float64) (a *Adaptive) {
	a = &Adaptive{
		np:        np,
		pnames:    pnames,
		t:         make([]int, np),
		loct:      make([]int, np),
		bmu:       make([]float64, np),
		mu:        make([]float64, np),
		vals:      make([]chan float64, np),
		sum:       make([]float64, np),
		sumsq:     make([]float64, np),
		SD:        make([]float64, np),
		delta:     make([]bool, np),
		converged: make([]bool, np),
	}
	a.K = 20
	a.Skip = 2000
	a.MaxUpdate = 2000
	a.Epsilon = 5e-1
	a.C = 1
	a.Nu = 3
	a.Lambda = 2.4

	for p := 0; p < np; p++ {
		a.mu[p] = math.NaN()
		a.vals[p] = make(chan float64, a.wSize)
		a.SD[p] = sd
	}

	return
}

func (a *Adaptive) SetWSize(wSize int) {
	a.wSize = wSize
	for p:=0 ; p < a.np; p ++ {
		a.vals[p] = make(chan float64, a.wSize)
	}
}

func (a *Adaptive) RobbinsMonro(p int) (gamma, tdelta float64) {
	tdelta = a.bmu[p] - a.mu[p]
	if (tdelta > 0 && !a.delta[p]) || (tdelta < 0 && a.delta[p]) {
		a.loct[p]++
	}
	a.delta[p] = tdelta > 0
	beta := 1 / math.Max(1, 1+a.Nu)
	gamma = a.C * math.Pow(float64(a.loct[p]), beta)
	return
}

func (a *Adaptive) CheckConvergenceMu(p int, val float64) {
	if len(a.vals[p]) == a.wSize {
		oldVal := <-a.vals[p]
		a.sum[p] -= oldVal
		a.sumsq[p] -= oldVal * oldVal
	}
	a.vals[p] <- val
	a.sum[p] += val
	a.sumsq[p] += val * val
	if len(a.vals[p]) == a.wSize {
		mean := a.sum[p] / float64(len(a.vals[p]))
		a.SD[p] = math.Sqrt(a.sumsq[p]/float64(len(a.vals[p])) - mean*mean)
		if a.SD[p]/mean < a.Epsilon || a.t[p]/a.K > a.MaxUpdate {
			a.converged[p] = true
			a.nconverged++
			var reason string
			if a.SD[p]/mean < a.Epsilon {
				reason = "SD/mean"
				log.Print("(reason: SD/mean)")
			} else {
				reason = "max_update"
			}
			log.Printf("%s converged, reason: %s (%d/%d)", a.pnames[p], reason, a.nconverged, a.np)
		}
	}
}

func (a *Adaptive) UpdateMu(p int, val float64) {
	if a.converged[p] {
		return
	}
	if math.IsNaN(a.mu[p]) {
		a.mu[p] = val
	}

	// Recursive mu formua
	a.bmu[p] = val/float64(a.K) + a.bmu[p]

	if a.t[p] > 0 && a.t[p]%a.K == 0 {
		gamma, delta := a.RobbinsMonro(p)

		// reset batch mu
		a.bmu[p] = 0
		a.mu[p] = a.mu[p] + gamma*delta
		a.CheckConvergenceMu(p, val)
	}
	a.t[p]++
}

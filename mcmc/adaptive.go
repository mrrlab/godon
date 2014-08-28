package mcmc

import (
	"fmt"
	"log"
	"math"
	"math/rand"
)

type AdaptiveParameter struct {
	*Float64Parameter
	t         int
	loct      int
	locMean   float64
	m2        float64
	bmu       float64
	mu        float64
	bvariance float64
	variance  float64
	vals      chan float64
	delta     bool
	converged bool

	*AdaptiveSettings
}

type AdaptiveSettings struct {
	WSize     int
	K         int
	Skip      int
	MaxUpdate int
	Epsilon   float64
	C         float64
	Nu        float64
	Lambda    float64
	SD        float64
}

func square(x float64) float64 {
	return x * x
}

func NewAdaptiveSettings() *AdaptiveSettings {
	return &AdaptiveSettings{
		WSize:     10,
		K:         20,
		Skip:      2000,
		MaxUpdate: 2000,
		Epsilon:   5e-1,
		C:         1,
		Nu:        3,
		Lambda:    2.4,
		SD:        1e-2,
	}
}

func NewAdaptiveParameter(par *float64, name string, ap *AdaptiveSettings) (a *AdaptiveParameter) {
	a = &AdaptiveParameter{
		Float64Parameter: NewFloat64Parameter(par, name),
		AdaptiveSettings: ap,
	}
	a.mu = math.NaN()
	a.vals = make(chan float64, a.WSize)
	if a.SD <= 0 {
		panic("SD should be >= 0")
	}
	a.variance = square(a.SD)

	a.ProposalFunc = a.AdaptiveProposal()

	return
}

func (a *AdaptiveSettings) String() string {
	return fmt.Sprintf("Adaptive MCMC settings <WSize=%v, K=%v, Skip=%v, MaxUpdate=%v, C=%v, Nu=%v, Lambda=%v>",
		a.WSize, a.K, a.Skip, a.MaxUpdate, a.C, a.Nu, a.Lambda)
}

func (a *AdaptiveParameter) Accept() {
	a.UpdateMu()
}

func (a *AdaptiveParameter) RobbinsMonro() (gamma, udelta, vdelta float64) {
	udelta = a.bmu - a.mu
	vdelta = a.bvariance - a.variance
	if (udelta > 0 && !a.delta) || (udelta < 0 && a.delta) {
		a.loct++
	}
	a.delta = udelta > 0
	beta := 1 / math.Max(1, 1+a.Nu)
	gamma = a.C * math.Pow(float64(a.loct), beta)
	return
}

func (a *AdaptiveParameter) CheckConvergenceMu() {
	if len(a.vals) == a.WSize {
		oldVal := <-a.vals
		delta := oldVal - a.locMean
		a.locMean -= delta / float64(len(a.vals))
		a.m2 -= delta * (oldVal - a.locMean)
	}

	a.vals <- *a.float64
	delta := *a.float64 - a.locMean
	a.locMean += delta / float64(len(a.vals))
	a.m2 += delta * (*a.float64 - a.locMean)

	variance := a.m2 / float64(len(a.vals)-1)

	if len(a.vals) == a.WSize {
		sd := math.Sqrt(variance)
		if sd/a.locMean < a.Epsilon || a.t/a.K > a.MaxUpdate {
			a.converged = true
			var reason string
			if sd/a.locMean < a.Epsilon {
				reason = "SD/mean"
				log.Print("(reason: SD/mean)")
			} else {
				reason = "max_update"
			}
			log.Printf("%s converged, reason: %s", a.Name(), reason)
		}
	}
}

func (a *AdaptiveParameter) UpdateMu() {
	if a.converged {
		return
	}
	if math.IsNaN(a.mu) {
		a.mu = *a.float64
	}

	// Recursive mu formua
	a.bmu = *a.float64/float64(a.K) + a.bmu
	a.bvariance = square(*a.float64)/float64(a.K-1) + a.bvariance

	if a.t > 0 && a.t%a.K == 0 {
		gamma, udelta, vdelta := a.RobbinsMonro()

		// reset batch mu
		a.bmu = 0
		a.mu += gamma * udelta
		a.variance += gamma * vdelta
		a.CheckConvergenceMu()
	}
	a.t++
}

func (a *AdaptiveParameter) AdaptiveProposal() func(float64) float64 {
	return func(x float64) float64 {
		return x + rand.NormFloat64()*math.Sqrt(a.variance)
	}
}

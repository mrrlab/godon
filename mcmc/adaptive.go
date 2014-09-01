package mcmc

import (
	"fmt"
	"log"
	"math"
	"math/rand"
)

type AdaptiveParameter struct {
	*Float64Parameter
	iter int
	t    int
	loct int

	//parameters
	mean     float64
	variance float64
	delta    bool

	//batch parameters
	bmean float64
	bm2   float64

	//convergence check
	vals      chan float64
	cmean     float64
	cm2       float64
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
		Skip:      500,
		MaxUpdate: 200,
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
	a.mean = math.NaN()
	a.vals = make(chan float64, a.WSize)
	if a.SD <= 0 {
		panic("SD should be >= 0")
	}
	if a.K < 2 {
		panic("K should be >= 2")
	}
	a.variance = square(a.SD)

	a.ProposalFunc = a.AdaptiveProposal()

	return
}

func (a *AdaptiveSettings) String() string {
	return fmt.Sprintf("Adaptive MCMC settings <WSize=%v, K=%v, Skip=%v, MaxUpdate=%v, C=%v, Nu=%v, Lambda=%v>",
		a.WSize, a.K, a.Skip, a.MaxUpdate, a.C, a.Nu, a.Lambda)
}

func (a *AdaptiveParameter) Accept(iter int) {
	if iter >= a.Skip {
		a.UpdateMu()
	}
}

func (a *AdaptiveParameter) RobbinsMonro() (gamma float64) {
	delta := a.bmean - a.mean
	if (delta > 0 && !a.delta) || (delta < 0 && a.delta) {
		a.loct++
	}
	a.delta = delta > 0
	beta := 1 / math.Max(1, 1+a.Nu)
	gamma = a.C / math.Pow(float64(a.loct+1), beta)
	return
}

func (a *AdaptiveParameter) CheckConvergenceMu() {
	if len(a.vals) == a.WSize {
		oldVal := <-a.vals
		delta := oldVal - a.cmean
		a.cmean -= delta / float64(len(a.vals))
		a.cm2 -= delta * (oldVal - a.cmean)
	}

	a.vals <- *a.float64
	delta := *a.float64 - a.cmean
	a.cmean += delta / float64(len(a.vals))
	a.cm2 += delta * (*a.float64 - a.cmean)

	if len(a.vals) == a.WSize {
		variance := a.cm2 / float64(len(a.vals)-1)
		sd := math.Sqrt(variance)
		if sd/a.cmean < a.Epsilon || a.t/a.K > a.MaxUpdate {
			a.converged = true
			var reason string
			if sd/a.cmean < a.Epsilon {
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
	if math.IsNaN(a.mean) {
		a.mean = *a.float64
	}
	// Incremental batch mean and variance
	// index in batch 0 .. a.K-1
	bi := (a.t - a.Skip) % a.K

	if (a.t-a.Skip) > 0 && bi == 0 {
		gamma := a.RobbinsMonro()

		bvariance := a.bm2 / float64(a.K-1)

		a.mean += gamma * (a.bmean - a.mean)
		a.variance += gamma * (bvariance - a.variance)

		a.CheckConvergenceMu()

		// reset batch mu
		a.bmean = 0
		a.bm2 = 0
	}

	delta := *a.float64 - a.bmean
	a.bmean += delta / float64(bi+1)
	a.bm2 += delta * (*a.float64 - a.bmean)
	// there is no need to calculate this every iterations

	a.t++
}

func (a *AdaptiveParameter) AdaptiveProposal() func(float64) float64 {
	return func(x float64) float64 {
		return x + rand.NormFloat64()*math.Sqrt(a.variance)*a.Lambda
	}
}

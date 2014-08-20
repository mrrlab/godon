package mcmc

import (
	"fmt"
	"log"
	"math"
	"math/rand"
)

const (
	K          = 20
	W_SIZE     = 10
	MAX_UPDATE = 200
	EPSILON2   = 5e-1
	C          = 1
	NU         = 3
	LAMBDA     = 2.4
)

func AMCMC(m Optimizable, burnIn, iterations int, report int) {
	accepted := 0
	np := m.GetNumberOfParameters()
	L := m.Likelihood()
	if burnIn > 0 {
		log.Printf("Burnin for %d iterations", burnIn)
	}
	fmt.Printf("iteration likelihood %s\n", ParameterNamesString(m))

	t := make([]int, np)
	loct := make([]int, np)
	bmu := make([]float64, np)
	mu := make([]float64, np)
	vals := make([]chan float64, np)
	sum := make([]float64, np)
	sumsq := make([]float64, np)
	stdev := make([]float64, np)
	for p, _ := range mu {
		mu[p] = m.GetParameter(p)
		vals[p] = make(chan float64, W_SIZE)
		stdev[p] = STDEV
	}
	delta := make([]bool, np)
	converged := make([]bool, np)
	nconverged := 0

	for i := 0; i < burnIn+iterations+burnIn; i++ {
		if i%ACCEPTED_RATE == 0 {
			log.Printf("Acceptance rate %f%%", 100*float64(accepted)/float64(ACCEPTED_RATE))
			accepted = 0
		}

		iter := i - burnIn
		if iter == 0 {
			log.Print("Starting sampling")
		}
		if iter >= 0 && iter%report == 0 {
			log.Printf("%d: L=%f", i, L)
			fmt.Printf("%d %f %s\n", i, L, ParameterString(m))
		}
		p := rand.Intn(np)
		val := m.GetParameter(p)
		newVal := val + rand.NormFloat64()*stdev[p]*LAMBDA
		m.SetParameter(p, newVal)
		newL := m.Likelihood()
		a := math.Exp(newL - L)
		if a < 1 && rand.Float64() > a {
			m.SetParameter(p, val)
		} else {
			if !converged[p] {
				bmu[p] = newVal/K + bmu[p]
				if t[p] > 0 && t[p] > 0 && t[p]%K == 0 {
					tdelta := bmu[p] - mu[p]
					if (tdelta > 0 && !delta[p]) || (tdelta < 0 && delta[p]) {
						loct[p]++
					}
					delta[p] = tdelta > 0
					beta := 1 / math.Max(1, 1+NU)
					gamma := C * math.Pow(float64(loct[p]), beta)
					bmu[p] = 0
					mu[p] = mu[p] + gamma*tdelta
					if len(vals[p]) == W_SIZE {
						oldVal := <-vals[p]
						sum[p] -= oldVal
						sumsq[p] -= oldVal * oldVal
					}
					vals[p] <- newVal
					sum[p] += newVal
					sumsq[p] += newVal * newVal
					if len(vals[p]) == W_SIZE {

						mean := sum[p] / float64(len(vals[p]))
						stdev[p] = math.Sqrt(sumsq[p]/float64(len(vals[p])) - mean*mean)
						log.Printf("%v stdev = %v", m.GetParameterName(p), stdev[p])
						if stdev[p]/mean < EPSILON2 || t[p]/K > MAX_UPDATE {
							converged[p] = true
							nconverged++
							log.Printf("%s converged (%d/%d)", m.GetParameterName(p), nconverged, np)
						}
					}
				}
				t[p]++
			}
			L = newL
			accepted++
		}
	}
	log.Print("Finished MCMC")
	for i := 0; i < np; i++ {
		log.Printf("%s=%f", m.GetParameterName(i), m.GetParameter(i))
	}
}

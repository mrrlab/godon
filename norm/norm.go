package main

import (
	"flag"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
	"syscall"
	"time"

	"bitbucket.org/Davydov/godon/optimize"
)

const (
	N    = 2
	K    = 100
	u    = 3
	sdev = 2
)

type MNormModel struct {
	data             [][]float64
	mean             []float64
	sd               []float64
	n                int
	adaptiveSettings *optimize.AdaptiveSettings
	parameters       optimize.FloatParameters
}

func square(x float64) float64 {
	return x * x
}

func NewMNormModel(data [][]float64, adaptive bool) (m *MNormModel) {
	mean := make([]float64, len(data))
	sd := make([]float64, len(data))
	parameters := make(optimize.FloatParameters, 0, len(data)*2)
	for i, _ := range sd {
		sd[i] = 1
	}
	m = &MNormModel{data: data,
		mean:       mean,
		sd:         sd,
		n:          len(data),
		parameters: parameters,
	}
	var nfp optimize.NewFloatParameter
	if adaptive {
		m.adaptiveSettings = optimize.NewAdaptiveSettings()
		nfp = m.adaptiveSettings.NewParameter
	} else {
		nfp = optimize.NewFloatParameterBasic
	}
	m.AddParameters(nfp)
	return
}

func (m *MNormModel) Copy() optimize.Optimizable {
	mean := make([]float64, m.n)
	copy(mean, m.mean)
	sd := make([]float64, m.n)
	copy(sd, m.sd)
	parameters := make(optimize.FloatParameters, 0, m.n*2)
	newM := &MNormModel{data: m.data,
		mean:             mean,
		sd:               sd,
		n:                m.n,
		parameters:       parameters,
		adaptiveSettings: m.adaptiveSettings,
	}
	var nfp optimize.NewFloatParameter
	if newM.adaptiveSettings != nil {
		nfp = m.adaptiveSettings.NewParameter
	} else {
		nfp = optimize.NewFloatParameterBasic
	}
	newM.AddParameters(nfp)
	return newM
}

func (m *MNormModel) AddParameters(nfp optimize.NewFloatParameter) {
	for i := 0; i < m.n; i++ {
		name := "sd" + strconv.Itoa(i)
		par := nfp(&m.sd[i], name)
		par.SetMin(0)
		par.SetMax(100)
		par.SetPriorFunc(optimize.UniformPrior(0, 100, false, false))
		par.SetProposalFunc(optimize.NormalProposal(0.1))
		m.parameters.Append(par)

		name = "mean" + strconv.Itoa(i)
		par = nfp(&m.mean[i], name)
		par.SetMin(-100)
		par.SetMax(100)
		par.SetPriorFunc(optimize.UniformPrior(-100, 100, false, false))
		par.SetProposalFunc(optimize.NormalProposal(0.1))
		m.parameters.Append(par)
	}
}

func (m *MNormModel) GetFloatParameters() optimize.FloatParameters {
	return m.parameters
}

func (m *MNormModel) Likelihood() (res float64) {
	for i := 0; i < m.n; i++ {
		if m.sd[i] < 0 {
			return math.Inf(-1)
		}
		for _, x := range m.data[i] {
			lnL := -math.Log(m.sd[i]*math.Sqrt(2*math.Pi)) - square(x-m.mean[i])/2/square(m.sd[i])
			res += lnL
		}
	}
	return
}

func genData(mean []float64, sd []float64, n int) (data [][]float64) {
	if len(mean) != len(sd) {
		panic("sdev and mean should be the same length")
	}
	npar := len(mean)

	data = make([][]float64, npar)

	for i := 0; i < npar; i++ {
		data[i] = make([]float64, 0, n)
		for j := 0; j < n; j++ {
			data[i] = append(data[i], rand.NormFloat64()*sd[i]+mean[i])
		}
	}

	return data

}

func getMeanSD(data []float64) (mean, sd float64) {
	for _, x := range data {
		mean += x
	}
	mean /= float64(len(data))
	for _, x := range data {
		sd += (mean - x) * (mean - x)
	}
	sd /= float64(len(data) - 1)
	sd = math.Sqrt(sd)
	return
}

func main() {
	amcmc := flag.Bool("a", false, "adaptive mcmc")
	simplex := flag.Bool("simplex", false, "use downhill simplex method")
	iter := flag.Int("iter", 100000, "number of iterations")
	seed := flag.Int64("seed", -1, "random generator seed, default time based")
	flag.Parse()

	log.Print("Starting")
	log.Printf("Will generate %d values for %d normal distributions", K, N)
	if *seed == -1 {
		rand.Seed(time.Now().UnixNano())
	} else {
		rand.Seed(*seed)
	}
	mean := make([]float64, 0, K)
	sd := make([]float64, 0, K)

	for i := 0; i < N; i++ {
		mean = append(mean, float64(rand.Intn(100)-50))
		sd = append(sd, float64(rand.Intn(10)+1))
	}

	data := genData(mean, sd, K)

	for i, xs := range data {
		m, s := getMeanSD(xs)
		log.Printf("%v: Norm(%v, %v^2), mean=%v, sd=%v", i, mean[i], sd[i], m, s)
	}

	m := NewMNormModel(data, *amcmc)

	var opt optimize.Optimizer
	if !*simplex {
		chain := optimize.NewMH(false, 0)
		chain.AccPeriod = 200
		opt = chain
	} else {
		simp := optimize.NewDS()
		opt = simp
	}

	opt.SetReportPeriod(100)
	opt.SetOptimizable(m)
	opt.WatchSignals(os.Interrupt, syscall.SIGUSR2)
	opt.Run(*iter)
}

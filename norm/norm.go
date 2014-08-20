package main

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"time"

	"bitbucket.com/Davydov/golh/mcmc"
)

const (
	N    = 2
	K    = 100
	u    = 3
	sdev = 2
)

type Optimizable interface {
	SetDefaults()
	GetNumberOfParameters() int
	GetParameterName(i int) string
	GetParameter(i int) float64
	SetParameter(i int, val float64)
	Likelihood() float64
}

type MNormModel struct {
	data [][]float64
	mean []float64
	sd   []float64
	n    int
}

func square(x float64) float64 {
	return x * x
}
func NewMNormModel(data [][]float64) *MNormModel {
	mean := make([]float64, len(data))
	sd := make([]float64, len(data))
	for i, _ := range sd {
		sd[i] = 1
	}
	return &MNormModel{data: data,
		mean: mean,
		sd:   sd,
		n:    len(data)}
}

func (m *MNormModel) SetDefaults() {
}

func (m *MNormModel) GetNumberOfParameters() int {
	return m.n * 2
}

func (m *MNormModel) GetParameterName(i int) (s string) {
	if i >= 2*m.n {
		panic("Unknown parameter")
	}
	if i < m.n {
		s = "mean"
	} else {
		s = "sd"
	}
	return fmt.Sprintf("%s%d", s, i%m.n)
}

func (m *MNormModel) GetParameter(i int) float64 {
	if i >= 2*m.n {
		panic("Unknown parameter")
	}
	if i < m.n {
		return m.mean[i]
	} else {
		return m.sd[i%m.n]
	}
}

func (m *MNormModel) SetParameter(i int, val float64) {
	if i >= 2*m.n {
		panic("Unknown parameter")
	}
	if i < m.n {
		m.mean[i] = val
	} else {
		val = math.Max(math.Abs(val), 0)
		m.sd[i%m.n] = val
	}
}

func (m *MNormModel) Likelihood() (res float64) {
	for i := 0; i < m.n; i++ {
		for _, x := range m.data[i] {
			lnL := - math.Log(m.sd[i] * math.Sqrt(2*math.Pi)) - square(x-m.mean[i])/2/square(m.sd[i])
			res += lnL
		}
	}
	return
}

func init() {
	rand.Seed(time.Now().UnixNano())
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
	log.Print("Starting")
	log.Printf("Will generate %d values for %d normal distributions", K, N)
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

	m := NewMNormModel(data)

	mcmc.MCMC(m, 0, 100000, 1000)
}

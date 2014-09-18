package optimize

import (
	"log"
	"math"
)

const (
	TINY        = 1e-10
	SMALL_DELTA = 1.1
)

type DS struct {
	BaseOptimizer
	delta         float64
	ftol          float64
	repeat        bool
	oldL          float64
	points        []Optimizable
	psum          []float64
	allparameters []FloatParameters
	l             []float64
}

func NewDS() (ds *DS) {
	ds = &DS{
		delta: 1,
		ftol:  1E-10,
	}
	return
}

func (ds *DS) createSimplex(delta float64) {
	ds.points = make([]Optimizable, len(ds.parameters)+1)
	ds.allparameters = make([]FloatParameters, len(ds.points))
	ds.l = make([]float64, len(ds.points))
	for i := range ds.points {
		point := ds.Optimizable.Copy()
		ds.points[i] = point
		ds.allparameters[i] = point.GetFloatParameters()
	}
	for i := 0; i < len(ds.parameters); i++ {
		parameter := ds.allparameters[i+1][i]
		parameter.Set(parameter.Get() + delta)
	}
	for i := range ds.points {
		if ds.allparameters[i].InRange() {
			ds.l[i] = ds.points[i].Likelihood()
		} else {
			ds.l[i] = math.Inf(-1)
			println("Warning: out of range", i)
		}
	}
}

// amotry extrapolates by factor fac throught the face of the simplex accros from
// the low point, tries it, and replaces the high point if the new point is better.
func (ds *DS) amotry(ilo int, fac float64) float64 {
	ds.calcPsum()
	ndim := len(ds.parameters)
	fac1 := (1 - fac) / float64(ndim)
	fac2 := fac1 - fac
	for j := 0; j < ndim; j++ {
		ds.parameters[j].Set(ds.psum[j]*fac1 - ds.allparameters[ilo][j].Get()*fac2)
	}
	var l float64
	if ds.parameters.InRange() {
		l = ds.Likelihood()
	} else {
		l = math.Inf(-1)
	}
	if l > ds.l[ilo] {
		ds.points[ilo], ds.Optimizable = ds.Optimizable, ds.points[ilo]
		ds.allparameters[ilo], ds.parameters = ds.parameters, ds.allparameters[ilo]
		ds.l[ilo] = l
	}
	return l
}

func (ds *DS) calcPsum() {
	ds.psum = make([]float64, len(ds.parameters))
	for i := range ds.psum {
		for _, parameters := range ds.allparameters {
			ds.psum[i] += parameters[i].Get()
		}
	}
}

func (ds *DS) Run(iterations int) {
	ds.createSimplex(ds.delta)
	// Lowest (worst), next-lowest and highest points
	var ilo, inlo, ihi int
	var llo, lnlo, lhi float64
	ds.maxL = math.Inf(-1)
Iter:
	for ds.i = 1; ds.i <= iterations; ds.i++ {
		if ds.l[0] < ds.l[1] {
			ilo = 0
			inlo = 1
			ihi = 1
		} else {
			ilo = 1
			inlo = 0
			ihi = 0
		}
		llo = ds.l[ilo]
		lnlo = ds.l[inlo]
		lhi = ds.l[ihi]
		for i := 2; i < len(ds.points); i++ {
			if ds.l[i] >= lhi {
				lhi = ds.l[i]
				ihi = i
			}
			if ds.l[i] < llo {
				lnlo = llo
				inlo = ilo
				llo = ds.l[i]
				ilo = i
			} else if ds.l[i] < lnlo {
				lnlo = ds.l[i]
				inlo = i
			}
		}
		if lhi > ds.maxL {
			ds.maxL = lhi
			ds.maxLPar = ds.allparameters[ihi].Values(ds.maxLPar)
		}
		ds.maxL = math.Max(ds.maxL, lhi)
		ds.BaseOptimizer.l = lhi
		_ = inlo
		if !ds.Quiet && ds.i%ds.repPeriod == 0 {
			log.Printf("%d: L=%f (%f)", ds.i, lhi, lhi-llo)
			ds.PrintLine(ds.allparameters[ihi], lhi)
			/*
				for i, parameters := range ds.allparameters {
					ds.PrintLine(parameters, ds.l[i])
				}
			*/
		}
		rtol := 2 * math.Abs(ds.l[ihi]-ds.l[ilo]) / (math.Abs(ds.l[ilo]) + math.Abs(ds.l[ihi]) + TINY)
		if rtol < ds.ftol {
			if ds.repeat && math.Abs(ds.oldL-lhi) < 2*TINY {
				ds.l[0], ds.l[ihi] = ds.l[ihi], ds.l[0]
				ds.points[0], ds.points[ihi] = ds.points[ihi], ds.points[0]
				ds.allparameters[0], ds.allparameters[ihi] = ds.allparameters[ihi], ds.allparameters[0]
				ds.Optimizable = ds.points[0].Copy()
				break Iter
			} else {
				ds.repeat = true
				ds.oldL = lhi
				ds.Optimizable = ds.points[ihi]
				ds.parameters = ds.allparameters[ihi]
				println(ds.parameters.ValuesString())
				ds.createSimplex(ds.delta)
				log.Printf("converged. retrying")
				continue
			}
		}
		l := ds.amotry(ilo, -1)
		switch {
		case l >= lhi:
			ds.amotry(ilo, 2)
		case l <= lnlo:
			lsave := llo
			l := ds.amotry(ilo, 0.5)
			if l <= lsave {
				for i, point := range ds.points {
					if i != ihi {
						for j := range ds.parameters {
							ds.allparameters[i][j].Set(0.5 * (ds.allparameters[i][j].Get() + ds.allparameters[ihi][j].Get()))
						}
						ds.l[i] = point.Likelihood()
					}

				}
			}
		}
		select {
		case s := <-ds.sig:
			log.Printf("Received signal %v, exiting.", s)
			break Iter
		default:
		}
	}
	if !ds.Quiet && ds.i == iterations {
		log.Printf("Iterations exceeded (%d)", iterations)
	}
	if !ds.Quiet {
		log.Print("Finished downhill simplex")
		log.Printf("Maximum likelihood: %v", ds.maxL)
		log.Printf("Parameter  names: %v", ds.parameters.NamesString())
		log.Printf("Parameter values: %v", ds.GetMaxLParameters())
		ds.PrintFinal()
	}

}

/*
func (ds *DS) Run() {

	m.L = m.Likelihood()
	m.MaxL = m.L
	m.MaxPar = m.ParameterString()
	m.PrintHeader()
	accepted := 0
Iter:
	for m.i = 0; m.i < iterations; m.i++ {
		if !m.Quiet && m.i > 0 && m.i%m.AccPeriod == 0 {
			log.Printf("Acceptance rate %.2f%%", 100*float64(accepted)/float64(m.AccPeriod))
			accepted = 0
		}

		if !m.Quiet && m.i%m.RepPeriod == 0 {
			log.Printf("%d: L=%f", m.i, m.L)
			m.PrintLine()
		}
		p := rand.Intn(len(m.parameters))
		par := m.parameters[p]
		par.Propose()
		newL := m.Likelihood()

		a := math.Exp(par.Prior() - par.OldPrior() + newL - m.L)
		if a > 1 || rand.Float64() < a {
			m.L = newL
			par.Accept(m.i)
			accepted++
			if m.L > m.MaxL {
				m.MaxL = m.L
				m.MaxPar = m.ParameterString()
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
		log.Printf("Maximum likelihood: %v", m.MaxL)
		log.Printf("Parameter  names: %v", m.ParameterNamesString())
		log.Printf("Parameter values: %v", m.MaxPar)
	}
	m.PrintFinal()
}

func (m *MH) PrintHeader() {
	if !m.Quiet {
		fmt.Printf("iteration\tlikelihood\t%s\n", m.ParameterNamesString())
	}
}

func (m *MH) PrintLine() {
	if !m.Quiet {
		fmt.Printf("%d\t%f\t%s\n", m.i, m.L, m.ParameterString())
	}
}

func (m *MH) PrintFinal() {
	if !m.Quiet {
		for _, par := range m.parameters {
			log.Printf("%s=%v", par.Name(), par.GetValue())
		}
	}
}

func (m *MH) ParameterNamesString() (s string) {
	for i, par := range m.parameters {
		if i != 0 {
			s += "\t"
		}
		s += par.Name()
	}
	return
}

func (m *MH) ParameterString() (s string) {
	for i, par := range m.parameters {
		if i != 0 {
			s += "\t"
		}
		s += par.GetValue().String()
	}
	return
}
*/

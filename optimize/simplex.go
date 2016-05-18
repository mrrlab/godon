package optimize

import (
	"math"
)

const (
	TINY        = 1e-10
	SMALL       = 1e-6
	SMALL_DELTA = 1.1
)

type DS struct {
	BaseOptimizer
	delta      float64
	ftol       float64
	repeat     bool
	oldL       float64
	points     []Optimizable
	psum       []float64
	parameters []FloatParameters
	l          []float64
	newOpt     Optimizable
	newPar     FloatParameters
}

func NewDS() (ds *DS) {
	ds = &DS{
		delta: 1,
		ftol:  TINY,
	}
	ds.repPeriod = 10
	return
}

func (ds *DS) createSimplex(opt Optimizable, delta float64) {
	parameters := opt.GetFloatParameters()
	ds.points = make([]Optimizable, len(parameters)+1)
	ds.parameters = make([]FloatParameters, len(ds.points))
	ds.l = make([]float64, len(ds.points))
	ds.points[0] = opt
	ds.parameters[0] = parameters
	for i := 1; i < len(ds.points); i++ {
		point := opt.Copy()
		ds.points[i] = point
		ds.parameters[i] = point.GetFloatParameters()
	}
	for i := 0; i < len(parameters); i++ {
		parameter := ds.parameters[i+1][i]
		parameter.Set(parameter.Get() + delta)
	}
	for i := range ds.points {
		if ds.parameters[i].InRange() {
			ds.l[i] = ds.points[i].Likelihood()
		} else {
			ds.l[i] = math.Inf(-1)
		}
	}
}

// amotry extrapolates by factor fac throught the face of the simplex accros from
// the low point, tries it, and replaces the high point if the new point is better.
func (ds *DS) amotry(ilo int, fac float64) float64 {
	if ds.newOpt == nil {
		ds.newOpt = ds.points[0].Copy()
		ds.newPar = ds.newOpt.GetFloatParameters()
	}
	ds.calcPsum()
	ndim := len(ds.newPar)
	fac1 := (1 - fac) / float64(ndim)
	fac2 := fac1 - fac
	for j := 0; j < ndim; j++ {
		ds.newPar[j].Set(ds.psum[j]*fac1 - ds.parameters[ilo][j].Get()*fac2)
	}
	var l float64
	if ds.newPar.InRange() {
		l = ds.newOpt.Likelihood()
	} else {
		l = math.Inf(-1)
	}
	if l > ds.l[ilo] {
		ds.points[ilo], ds.newOpt = ds.newOpt, ds.points[ilo]
		ds.parameters[ilo], ds.newPar = ds.newPar, ds.parameters[ilo]
		ds.l[ilo] = l
	}
	return l
}

func (ds *DS) calcPsum() {
	ds.psum = make([]float64, len(ds.parameters[0]))
	for i := range ds.psum {
		for _, parameters := range ds.parameters {
			ds.psum[i] += parameters[i].Get()
		}
	}
}

func (ds *DS) SetOptimizable(opt Optimizable) {
	ds.createSimplex(opt, ds.delta)
}

func (ds *DS) Run(iterations int) {
	// Lowest (worst), next-lowest and highest points
	var ilo, inlo, ihi int
	var llo, lnlo, lhi float64
	ds.PrintHeader(ds.parameters[0])
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
			ds.maxLPar = ds.parameters[ihi].Values(ds.maxLPar)
		}
		ds.maxL = math.Max(ds.maxL, lhi)
		ds.BaseOptimizer.l = lhi
		_ = inlo
		if ds.i%ds.repPeriod == 0 {
			log.Debugf("%d: L=%f (%f)", ds.i, lhi, lhi-llo)
			ds.PrintLine(ds.parameters[ihi], lhi)
			/*
				for i, parameters := range ds.allparameters {
					ds.PrintLine(parameters, ds.l[i])
				}
			*/
		}
		rtol := 2 * math.Abs(ds.l[ihi]-ds.l[ilo]) / (math.Abs(ds.l[ilo]) + math.Abs(ds.l[ihi]) + TINY)
		if rtol < ds.ftol {
			if ds.repeat && math.Abs(ds.oldL-lhi) < SMALL {
				break Iter
			} else {
				ds.repeat = true
				ds.oldL = lhi
				log.Infof("converged. retrying")
				ds.createSimplex(ds.points[ihi], ds.delta)
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
						for j := range ds.parameters[i] {
							ds.parameters[i][j].Set(0.5 * (ds.parameters[i][j].Get() + ds.parameters[ihi][j].Get()))
						}
						if ds.parameters[i].InRange() {
							ds.l[i] = point.Likelihood()
						} else {
							ds.l[i] = math.Inf(-1)
						}
					}

				}
			}
		}
		select {
		case s := <-ds.sig:
			log.Warningf("Received signal %v, exiting.", s)
			break Iter
		default:
		}
	}
	if ds.i == iterations {
		log.Warningf("Iterations exceeded (%d)", iterations)
	}

	log.Info("Finished downhill simplex")
	log.Noticef("Maximum likelihood: %v", lhi)
	log.Infof("Parameter  names: %v", ds.parameters[ihi].NamesString())
	log.Infof("Parameter values: %v", ds.parameters[ihi].ValuesString())
	ds.PrintFinal(ds.parameters[ihi])

}

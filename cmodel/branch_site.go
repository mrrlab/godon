package cmodel

import (
	"math"
	"math/big"
	"math/rand"
	"runtime"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/tree"
)

// BranchSite is an implementation of the branch-site model.
type BranchSite struct {
	*BaseModel
	q0, q1, q2             *codon.EMatrix
	kappa                  float64
	omega0, omega2         float64
	p01sum, p0prop         float64
	fixw2                  bool
	q0done, q1done, q2done bool
	propdone               bool
	summary                brachSiteSummary
}

// brachSiteSummary stores summary information.
type brachSiteSummary struct {
	SitePosteriorNEB []float64 `json:"sitePosteriorNEB,omitempty"`
	SitePosteriorBEB []float64 `json:"sitePosteriorBEB,omitempty"`
}

// NewBranchSite creates a new BranchSite model.
func NewBranchSite(cali codon.Sequences, t *tree.Tree, cf codon.Frequency, fixw2 bool) (m *BranchSite) {
	m = &BranchSite{
		fixw2: fixw2,
		q0:    &codon.EMatrix{CF: cf},
		q1:    &codon.EMatrix{CF: cf},
		q2:    &codon.EMatrix{CF: cf},
	}

	m.BaseModel = NewBaseModel(cali, t, cf, m)

	m.setupParameters()
	m.setBranchMatrices()
	m.SetDefaults()

	return

}

// GetNClass returns number of site classes.
func (m *BranchSite) GetNClass() int {
	return 4
}

// Copy makes a copy of the model preserving the model parameter
// values.
func (m *BranchSite) Copy() optimize.Optimizable {
	newM := &BranchSite{
		BaseModel: m.BaseModel.Copy(),
		q0:        &codon.EMatrix{CF: m.cf},
		q1:        &codon.EMatrix{CF: m.cf},
		q2:        &codon.EMatrix{CF: m.cf},
		kappa:     m.kappa,
		omega0:    m.omega0,
		omega2:    m.omega2,
		p01sum:    m.p01sum,
		p0prop:    m.p0prop,
		fixw2:     m.fixw2,
	}
	newM.BaseModel.Model = newM
	newM.setupParameters()
	newM.setBranchMatrices()
	return newM
}

// addParameters adds all the model parameters to the parameter
// storage.
func (m *BranchSite) addParameters(fpg optimize.FloatParameterGenerator) {
	kappa := fpg(&m.kappa, "kappa")
	kappa.SetOnChange(func() {
		m.q0done = false
		m.q1done = false
		m.q2done = false
	})
	kappa.SetPriorFunc(optimize.UniformPrior(0, 20, false, true))
	kappa.SetProposalFunc(optimize.NormalProposal(0.01))
	kappa.SetMin(1e-2)
	kappa.SetMax(20)
	m.parameters.Append(kappa)

	omega0 := fpg(&m.omega0, "omega0")
	omega0.SetOnChange(func() {
		m.q0done = false
	})
	omega0.SetPriorFunc(optimize.GammaPrior(1, 2, false))
	omega0.SetProposalFunc(optimize.NormalProposal(0.01))
	omega0.SetMin(1e-4)
	omega0.SetMax(1)
	m.parameters.Append(omega0)

	if !m.fixw2 {
		omega2 := fpg(&m.omega2, "omega2")
		omega2.SetOnChange(func() {
			m.q2done = false
		})
		omega2.SetPriorFunc(optimize.GammaPrior(1, 2, false))
		omega2.SetProposalFunc(optimize.NormalProposal(0.01))
		omega2.SetMin(1)
		omega2.SetMax(1000)
		m.parameters.Append(omega2)
	}

	p01sum := fpg(&m.p01sum, "p01sum")
	p01sum.SetOnChange(func() {
		m.propdone = false
	})
	p01sum.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
	p01sum.SetMin(0)
	p01sum.SetMax(1)
	p01sum.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(p01sum)

	p0prop := fpg(&m.p0prop, "p0prop")
	p0prop.SetOnChange(func() {
		m.propdone = false
	})
	p0prop.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
	p0prop.SetMin(0)
	p0prop.SetMax(1)
	p0prop.SetProposalFunc(optimize.NormalProposal(0.01))
	m.parameters.Append(p0prop)
}

// SetParameters sets the model parameter values.
func (m *BranchSite) SetParameters(kappa float64, omega0, omega2 float64, p0, p1 float64) {
	m.kappa = kappa
	m.omega0 = omega0
	if !m.fixw2 {
		m.omega2 = omega2
	} else {
		m.omega2 = 1
	}
	m.p01sum = p0 + p1
	m.p0prop = p0 / (p0 + p1)
	m.q0done, m.q1done, m.q2done = false, false, false
}

// GetParameters returns the model parameter values.
func (m *BranchSite) GetParameters() (kappa float64, omega0, omega2 float64, p0, p1 float64) {
	return m.kappa, m.omega0, m.omega2, m.p01sum * m.p0prop, m.p01sum * (1 - m.p0prop)
}

// SetDefaults sets the default initial parameter values.
func (m *BranchSite) SetDefaults() {
	// these parameters mostly come from codeml
	kappa := 1e-2 + rand.Float64()*10
	omega0 := 0.2 + 0.1*rand.Float64()
	omega2 := 3.1 + rand.Float64()
	x0 := 1.0 + 0.5*rand.Float64()
	x1 := 0.2 * rand.Float64()
	p0 := math.Exp(x0) / (1 + math.Exp(x0) + math.Exp(x1))
	p1 := math.Exp(x1) / (1 + math.Exp(x0) + math.Exp(x1))
	m.SetParameters(kappa, omega0, omega2, p0, p1)
}

// setBranchMatrices set matrices for all the branches.
func (m *BranchSite) setBranchMatrices() {
	for i := 0; i < len(m.qs); i++ {
		for _, node := range m.tree.NodeIdArray() {
			if node == nil {
				continue
			}
			switch i {
			case 0:
				m.qs[i][node.Id] = m.q0
			case 1:
				m.qs[i][node.Id] = m.q1
			case 2:
				if node.Class == 1 {
					m.qs[i][node.Id] = m.q2
				} else {
					m.qs[i][node.Id] = m.q0
				}
			case 3:
				if node.Class == 1 {
					m.qs[i][node.Id] = m.q2
				} else {
					m.qs[i][node.Id] = m.q1
				}
			}
		}
	}
}

// updateProportions updates proportions if model parameters are
// changing.
func (m *BranchSite) updateProportions() {
	p0 := m.p0prop * m.p01sum
	p1 := m.p01sum - p0
	m.prop[0][0] = p0
	m.prop[0][1] = p1
	m.prop[0][2] = (1 - p0 - p1) * p0 / (p0 + p1)
	m.prop[0][3] = (1 - p0 - p1) * p1 / (p0 + p1)

	for _, node := range m.tree.NodeIdArray() {
		if node == nil {
			continue
		}
		if node.Class == 1 {
			m.scale[node.Id] = m.prop[0][0]*m.q0.Scale + m.prop[0][1]*m.q1.Scale + (m.prop[0][2]+m.prop[0][3])*m.q2.Scale
		} else {
			m.scale[node.Id] = (m.prop[0][0]+m.prop[0][2])*m.q0.Scale + (m.prop[0][1]+m.prop[0][3])*m.q1.Scale
		}
	}
	m.propdone = true
	m.expAllBr = false
}

// updateMatrices updates matrices if model parameters are changing.
func (m *BranchSite) updateMatrices() {
	if !m.q0done {
		Q0, s0 := codon.CreateTransitionMatrix(m.cf, m.kappa, m.omega0, m.q0.Q)
		m.q0.Set(Q0, s0)
		err := m.q0.Eigen()
		if err != nil {
			panic("error eigen q0")
		}
		m.q0done = true
	}

	if !m.q1done {
		Q1, s1 := codon.CreateTransitionMatrix(m.cf, m.kappa, 1, m.q1.Q)
		m.q1.Set(Q1, s1)
		err := m.q1.Eigen()
		if err != nil {
			panic("error eigen q1")
		}
		m.q1done = true
	}

	if !m.q2done {
		Q2, s2 := codon.CreateTransitionMatrix(m.cf, m.kappa, m.omega2, m.q2.Q)
		m.q2.Set(Q2, s2)
		err := m.q2.Eigen()
		if err != nil {
			panic("error eigen q2")
		}
		m.q2done = true
	}

	m.propdone = false
}

//siteLMatrix returns site likelihood matrix. First index is w0, second
//w2, third class, forth position.
func (m *BranchSite) siteLMatrix(w0, w2 []float64) (res [][][][]float64) {
	res = make([][][][]float64, len(w0))
	nClass := m.GetNClass()
	nPos := m.cali.Length()
	nni := m.tree.MaxNodeId() + 1

	// temporary storage for likelihood computation
	plh := make([][]float64, nni)
	for i := 0; i < nni; i++ {
		plh[i] = make([]float64, m.cf.GCode.NCodon+1)
	}

	counter := 0

	nWorkers := runtime.GOMAXPROCS(0)
	done := make(chan struct{}, nWorkers)
	type bebtask struct {
		iW0, iW2, class, pos int
	}
	tasks := make(chan bebtask, nPos)

	go func() {
		nni := m.tree.MaxNodeId() + 1
		plh := make([][]float64, nni)
		for i := 0; i < nni; i++ {
			plh[i] = make([]float64, m.cf.GCode.NCodon+1)
		}
		for task := range tasks {
			res[task.iW0][task.iW2][task.class][task.pos] = m.fullSubL(task.class, task.pos, plh)
			done <- struct{}{}
		}
	}()

	for iW0, w0 := range w0 {
		m.omega0 = w0
		m.q0done = false
		res[iW0] = make([][][]float64, len(w2))
		for iW2, w2 := range w2 {
			m.omega2 = w2
			m.q2done = false
			m.updateMatrices()
			// in this scenario we keep q-factor as computed from MLE
			m.ExpBranches()
			res[iW0][iW2] = make([][]float64, nClass)
			for class := 0; class < nClass; class++ {
				switch {
				case class == 0 && iW2 != 0:
					res[iW0][iW2][class] = res[iW0][0][class]
				case class == 1 && (iW0 != 0 || iW2 != 0):
					res[iW0][iW2][class] = res[0][0][class]
				case class == 3 && iW0 != 0:
					res[iW0][iW2][class] = res[0][iW2][class]
				default:
					res[iW0][iW2][class] = make([]float64, nPos)

					counter++
					for pos := 0; pos < nPos; pos++ {
						//res[i_w0][i_w2][class][pos] = m.fullSubL(class, pos, plh)
						tasks <- bebtask{iW0, iW2, class, pos}
					}
					// wait for everyone to finish
					for pos := 0; pos < nPos; pos++ {
						<-done
					}
				}
			}
		}
	}
	log.Infof("Computed f(x_h|w) for %d classes", counter)
	return
}

// computePropBEB computes class proportions given the i, j and d parameters.
func (m *BranchSite) computePropBEB(prop []float64, i, j, d int) []float64 {
	if prop == nil {
		nClass := m.GetNClass()
		prop = make([]float64, nClass)
	}
	prop[0] = (1. + float64(j)/2*3 + float64(j%2)) / (3 * float64(d))
	prop[1] = (1. + float64(d-1-i)*3 + float64(j%2)) / (3 * float64(d))
	prop[2] = (1 - prop[0] - prop[1]) * prop[0] / (prop[0] + prop[1])
	prop[3] = (1 - prop[0] - prop[1]) * prop[1] / (prop[0] + prop[1])
	return prop
}

// BEBPosterior returns BEB posterior values.
func (m *BranchSite) BEBPosterior() (res []float64) {
	nPos := m.cali.Length()
	nClass := m.GetNClass()
	res = make([]float64, nPos)

	log.Info("w0 and w2 grid")

	// first create a grid of parameters
	w0s := floatRange(0.05, 0.1, 10)
	log.Infof("w0: %v", strFltSlice(w0s))

	w2s := floatRange(1.5, 1, 10)
	log.Infof("w2: %v", strFltSlice(w2s))

	matr := m.siteLMatrix(w0s, w2s)

	// normalization
	fx := big.NewFloat(0)

	// storing sum values for every position
	su := make([]float64, nPos)

	var prop []float64

	product := new(big.Float)
	productNoPos := new(big.Float)
	tmp := new(big.Float)

	posterior := make([]big.Float, nPos)

	d := 10
	// first compute the stat sum (fx)
	for iW0 := range w0s {
		for iW2 := range w2s {
			for i := 0; i < d; i++ {
				for j := 0; j <= 2*i; j++ {
					prop = m.computePropBEB(prop, i, j, d)
					// prior for p0 and p1 0.01, for w0 and w2 0.1
					prior := 0.01 * 0.1 * 0.1

					product.SetFloat64(1)
					// first compute the product for every site
					for pos := 0; pos < nPos; pos++ {
						su[pos] = 0
						for class := 0; class < nClass; class++ {
							su[pos] += matr[iW0][iW2][class][pos] * prop[class]
						}
						tmp.SetFloat64(su[pos])
						product.Mul(product, tmp)
					}
					tmp.SetFloat64(prior)
					product.Mul(product, tmp)

					fx.Add(fx, product)

					for pos := 0; pos < nPos; pos++ {
						// compute product of everything except pos
						tmp.SetFloat64(su[pos])
						productNoPos.Quo(product, tmp)

						tmp.SetFloat64(
							matr[iW0][iW2][2][pos]*prop[2] +
								matr[iW0][iW2][3][pos]*prop[3])
						tmp.Mul(tmp, productNoPos)
						posterior[pos].Add(&posterior[pos], tmp)
					}
				}
			}
		}
	}

	for pos := range posterior {
		posterior[pos].Quo(&posterior[pos], fx)
		res[pos], _ = posterior[pos].Float64()
	}

	return
}

// Final prints NEB results (only if with positive selection).
func (m *BranchSite) Final() {
	// if w2=1, do not perform NEB analysis.
	if m.fixw2 {
		log.Info("No NEB since no positive selection in the model.")
		return
	}
	classes := make(map[int]bool, m.GetNClass())
	classes[2] = true
	classes[3] = true

	posterior := m.NEBPosterior(classes)
	m.summary.SitePosteriorNEB = posterior

	log.Notice("NEB analysis")
	m.PrintPosterior(posterior)

	posterior = m.BEBPosterior()
	m.summary.SitePosteriorBEB = posterior

	log.Notice("BEB analysis")
	m.PrintPosterior(posterior)
}

// Likelihood computes likelihood.
func (m *BranchSite) Likelihood() float64 {
	if !m.q0done || !m.q1done || !m.q2done {
		m.updateMatrices()
	}
	if !m.propdone {
		m.updateProportions()
	}
	return m.BaseModel.Likelihood()
}

// Summary returns the run summary (site posterior for NEB and BEB).
func (m *BranchSite) Summary() interface{} {
	if m.summary.SitePosteriorBEB != nil || m.summary.SitePosteriorNEB != nil {
		return m.summary
	}
	// nil prevents json from printing "{}"
	return nil
}

package cmodel

import (
	"math"
	"math/big"
	"math/rand"
	"runtime"
	"time"

	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/dist"
	"bitbucket.org/Davydov/godon/optimize"
)

// BranchSiteGamma is an implementation of the branch-site model with gamma
// rates variation.
type BranchSiteGamma struct {
	*BaseModel
	q0s            []*codon.EMatrix
	q1s            []*codon.EMatrix
	q2s            []*codon.EMatrix
	kappa          float64
	omega0, omega2 float64
	p01sum, p0prop float64
	// site gamma alpha parameter
	alphas     float64
	gammas     []float64
	gammasprop []float64
	// parametrization inspired by Scheffler 2006
	// proportions
	ps1s, ps2s float64
	// rates, rs2s is inverted
	rs1s, rs2s float64
	// codon gamma alpha parameter
	alphac     float64
	gammac     []float64
	gammacprop []float64
	// parametrization inspired by Scheffler 2006
	// proportions
	ps1c, ps2c float64
	// rates, rs2c is inverted
	rs1c, rs2c float64
	// temporary array
	tmp                    []float64
	ncatcg                 int
	ncatsg                 int
	fixw2                  bool
	proportional           bool
	q0done, q1done, q2done bool
	propdone               bool
	gammacdone             bool
	gammasdone             bool
	summary                branchSiteGammaSummary
}

// branchSiteGammaSummary stores summary information.
type branchSiteGammaSummary struct {
	SitePosteriorNEB []float64 `json:"sitePosteriorNEB,omitempty"`
	SitePosteriorBEB []float64 `json:"sitePosteriorBEB,omitempty"`
	CodonGammaRates  []float64 `json:"codonGammaRates,omitempty"`
	PosteriorTime    float64   `json:"posteriorTime,omitempty"`
}

// Empty returns true if there's no data in the structure.
func (s branchSiteGammaSummary) Empty() bool {
	if s.SitePosteriorNEB == nil && s.SitePosteriorBEB == nil && s.CodonGammaRates == nil {
		return true
	}
	return false
}

// NewBranchSiteGamma creates a new BranchSiteGamma model.
func NewBranchSiteGamma(data *Data, fixw2 bool, ncatsg, ncatcg int, proportional bool) (m *BranchSiteGamma) {
	scat := ncatsg * ncatsg * ncatsg

	if proportional {
		if (ncatsg != 3 && ncatsg != 1) || (ncatcg != 3 && ncatcg != 1) {
			log.Fatal("Scheffler's parametrization is only supported for three discrete rates")
		}
	}
	m = &BranchSiteGamma{
		fixw2:        fixw2,
		proportional: proportional,
		ncatsg:       ncatsg,
		ncatcg:       ncatcg,
		q0s:          make([]*codon.EMatrix, scat*ncatcg),
		q1s:          make([]*codon.EMatrix, scat*ncatcg),
		q2s:          make([]*codon.EMatrix, scat*ncatcg),
		gammas:       make([]float64, ncatsg),
		gammasprop:   make([]float64, ncatsg),
		gammac:       make([]float64, ncatcg),
		gammacprop:   make([]float64, ncatcg),
		tmp:          make([]float64, maxInt(ncatcg, ncatsg, 3)),
	}
	m.BaseModel = NewBaseModel(data, m)

	for i := 0; i < scat*ncatcg; i++ {
		m.q0s[i] = codon.NewEMatrix(data.cFreq)
		m.q1s[i] = codon.NewEMatrix(data.cFreq)
		m.q2s[i] = codon.NewEMatrix(data.cFreq)
	}

	if !m.proportional {
		for i := range m.gammacprop {
			m.gammacprop[i] = 1 / float64(len(m.gammacprop))
		}
		for i := range m.gammasprop {
			m.gammasprop[i] = 1 / float64(len(m.gammasprop))
		}
	}

	m.setupParameters()
	m.setBranchMatrices()
	m.SetDefaults()

	return

}

// GetNClass returns number of site classes.
func (m *BranchSiteGamma) GetNClass() int {
	scat := m.ncatsg * m.ncatsg * m.ncatsg

	return 4 * scat * m.ncatcg
}

// Copy makes a copy of the model preserving the model parameter
// values.
func (m *BranchSiteGamma) Copy() optimize.Optimizable {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	newM := &BranchSiteGamma{
		BaseModel: m.BaseModel.Copy(),
		q0s:       make([]*codon.EMatrix, scat*m.ncatcg),
		q1s:       make([]*codon.EMatrix, scat*m.ncatcg),
		q2s:       make([]*codon.EMatrix, scat*m.ncatcg),
		kappa:     m.kappa,
		omega0:    m.omega0,
		omega2:    m.omega2,
		p01sum:    m.p01sum,
		p0prop:    m.p0prop,
		fixw2:     m.fixw2,

		proportional: m.proportional,

		ncatsg:     m.ncatsg,
		alphas:     m.alphas,
		gammas:     make([]float64, m.ncatsg),
		gammasprop: make([]float64, m.ncatsg),
		rs1s:       m.rs1s,
		rs2s:       m.rs2s,
		ps1s:       m.ps1s,
		ps2s:       m.ps2s,

		ncatcg:     m.ncatcg,
		alphac:     m.alphac,
		gammac:     make([]float64, m.ncatcg),
		gammacprop: make([]float64, m.ncatcg),
		rs1c:       m.rs1c,
		rs2c:       m.rs2c,
		ps1c:       m.ps1c,
		ps2c:       m.ps2c,

		tmp: make([]float64, maxInt(m.ncatcg, m.ncatsg, 3)),
	}
	newM.BaseModel.model = newM

	for i := 0; i < scat*m.ncatcg; i++ {
		newM.q0s[i] = codon.NewEMatrix(m.data.cFreq)
		newM.q1s[i] = codon.NewEMatrix(m.data.cFreq)
		newM.q2s[i] = codon.NewEMatrix(m.data.cFreq)
	}

	if !m.proportional {
		for i := range m.gammacprop {
			m.gammacprop[i] = 1 / float64(len(m.gammacprop))
		}
		for i := range m.gammasprop {
			m.gammasprop[i] = 1 / float64(len(m.gammasprop))
		}
	}

	newM.setupParameters()
	newM.setBranchMatrices()
	return newM
}

// addParameters adds all the model parameters to the parameter
// storage.
func (m *BranchSiteGamma) addParameters(fpg optimize.FloatParameterGenerator) {
	kappa := fpg(&m.kappa, "kappa")
	kappa.SetOnChange(func() {
		m.q0done = false
		m.q1done = false
		m.q2done = false
	})
	kappa.SetPriorFunc(optimize.UniformPrior(0, 20, false, true))
	kappa.SetProposalFunc(optimize.NormalProposal(0.01))
	kappa.SetMin(1e-2)
	kappa.SetMax(100)
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
	p01sum.SetMin(1e-12)
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

	if m.ncatsg > 1 {
		if !m.proportional {
			alphas := fpg(&m.alphas, "alphas")
			alphas.SetOnChange(func() {
				m.gammasdone = false
			})
			alphas.SetPriorFunc(optimize.GammaPrior(1, 2, false))
			alphas.SetMin(1e-2)
			alphas.SetMax(1000)
			alphas.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(alphas)
		} else {
			ps1s := fpg(&m.ps1s, "ps1s")
			ps1s.SetOnChange(func() {
				m.gammasdone = false
			})
			ps1s.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			ps1s.SetMin(1e-5)
			ps1s.SetMax(1 - 1e-5)
			ps1s.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(ps1s)
			ps2s := fpg(&m.ps2s, "ps2s")
			ps2s.SetOnChange(func() {
				m.gammasdone = false
			})
			ps2s.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			ps2s.SetMin(1e-5)
			ps2s.SetMax(1 - 1e-5)
			ps2s.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(ps2s)
			rs1s := fpg(&m.rs1s, "log_rs1s")
			rs1s.SetOnChange(func() {
				m.gammasdone = false
			})
			rs1s.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			rs1s.SetMin(-10)
			rs1s.SetMax(-1e-4)
			rs1s.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(rs1s)
			rs2s := fpg(&m.rs2s, "log_rs2s")
			rs2s.SetOnChange(func() {
				m.gammasdone = false
			})
			rs2s.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			rs2s.SetMin(1e-4)
			rs2s.SetMax(10)
			rs2s.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(rs2s)
		}
	}

	if m.ncatcg > 1 {
		if !m.proportional {
			alphac := fpg(&m.alphac, "alphac")
			alphac.SetOnChange(func() {
				m.gammacdone = false
			})
			alphac.SetPriorFunc(optimize.GammaPrior(1, 2, false))
			alphac.SetMin(1e-2)
			alphac.SetMax(1000)
			alphac.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(alphac)
		} else {
			ps1c := fpg(&m.ps1c, "ps1c")
			ps1c.SetOnChange(func() {
				m.gammacdone = false
			})
			ps1c.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			ps1c.SetMin(1e-5)
			ps1c.SetMax(1 - 1e-5)
			ps1c.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(ps1c)
			ps2c := fpg(&m.ps2c, "ps2c")
			ps2c.SetOnChange(func() {
				m.gammacdone = false
			})
			ps2c.SetPriorFunc(optimize.UniformPrior(0, 1, false, false))
			ps2c.SetMin(1e-5)
			ps2c.SetMax(1 - 1e-5)
			ps2c.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(ps2c)
			rs1c := fpg(&m.rs1c, "log_rs1c")
			rs1c.SetOnChange(func() {
				m.gammacdone = false
			})
			rs1c.SetPriorFunc(optimize.UniformPrior(-10, -1e-4, false, false))
			rs1c.SetMin(-10)
			rs1c.SetMax(-1e-4)
			rs1c.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(rs1c)
			rs2c := fpg(&m.rs2c, "log_rs2c")
			rs2c.SetOnChange(func() {
				m.gammacdone = false
			})
			rs2c.SetPriorFunc(optimize.UniformPrior(1e-4, 10, false, false))
			rs2c.SetMin(1e-4)
			rs2c.SetMax(10)
			rs2c.SetProposalFunc(optimize.NormalProposal(0.01))
			m.parameters.Append(rs2c)
		}
	}
}

// SetParameters sets the model parameter values.
func (m *BranchSiteGamma) SetParameters(kappa float64, omega0, omega2 float64, p0, p1 float64, alphas, alphac float64,
	rs1s, rs2s, ps1s, ps2s float64, rs1c, rs2c, ps1c, ps2c float64) {
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
	m.alphas = alphas
	m.alphac = alphac
	m.rs1s = rs1s
	m.rs2s = rs2s
	m.ps1s = ps1s
	m.ps2s = ps2s
	m.rs1c = rs1c
	m.rs2c = rs2c
	m.ps1c = ps1c
	m.ps2c = ps2c
}

// GetParameters returns the model parameter values.
func (m *BranchSiteGamma) GetParameters() (kappa float64, omega0, omega2 float64, p0, p1 float64, alphas, alphac float64,
	rs1s, rs2s, ps1s, ps2s float64, rs1c, rs2c, ps1c, ps2c float64) {
	return m.kappa, m.omega0, m.omega2, m.p01sum * m.p0prop, m.p01sum * (1 - m.p0prop), m.alphas, m.alphac, rs1s, rs2s, ps1s, ps2s, rs1c, rs2c, ps1c, ps2c
}

// SetDefaults sets the default initial parameter values.
func (m *BranchSiteGamma) SetDefaults() {
	// these parameters mostly come from codeml
	kappa := 1e-2 + rand.Float64()*10
	omega0 := 0.2 + 0.1*rand.Float64()
	omega2 := 3.1 + rand.Float64()
	x0 := 1.0 + 0.5*rand.Float64()
	x1 := 0.2 * rand.Float64()
	p0 := math.Exp(x0) / (1 + math.Exp(x0) + math.Exp(x1))
	p1 := math.Exp(x1) / (1 + math.Exp(x0) + math.Exp(x1))
	alphas := 0.5 + rand.Float64()*3
	alphac := 0.5 + rand.Float64()*3
	rs1s := -1e-4 - rand.Float64()*(10-1e-4)
	rs2s := 1e-4 + rand.Float64()*(10-1e-4)
	ps1s := 1e-5 + rand.Float64()*(1-2e-5)
	ps2s := 1e-5 + rand.Float64()*(1-2e-5)
	rs1c := -1e-4 - rand.Float64()*(10-1e-4)
	rs2c := 1e-4 + rand.Float64()*(10-1e-4)
	ps1c := 1e-5 + rand.Float64()*(1-2e-5)
	ps2c := 1e-5 + rand.Float64()*(1-2e-5)

	m.SetParameters(kappa, omega0, omega2, p0, p1, alphas, alphac,
		rs1s, rs2s, ps1s, ps2s, rs1c, rs2c, ps1c, ps2c)
}

// Organization of the class categories.
// cat 1, internal gamma cat 1, external gamma cat 1
// cat 1, internal gamma cat 1, external gamma cat 2
// ...
// cat 1, internal gamma cat 1, external gamma ncatcg
// ...
// cat 1, internal gamma cat ncatsg^3, external gamma cat 1
// ...
// cat 1, internal gamma cat ncatsg^3, external gamma ncatcg
// cat 2, internal gamma cat 1, external gamma cat 1
// ...
// cat 2, internal gamma cat ncatsg^3, external gamma cat ncatcg
// ...
// ...
// cat 4, internal gamma cat 1, external gamma cat 1
// ...
// cat 4, internal gamma cat ncatsg^3, external gamma cat ncatcg
// (total: 4 * ncatsg^3 * ncatcg)

// setBranchMatrices set matrices for all the branches.
func (m *BranchSiteGamma) setBranchMatrices() {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	for _, node := range m.data.Tree.NodeIDArray() {
		if node == nil {
			continue
		}
		for i := 0; i < m.ncatcg; i++ {
			for j := 0; j < scat; j++ {
				m.qs[i+j*m.ncatcg+0*scat*m.ncatcg][node.ID] = m.q0s[i+j*m.ncatcg]
				m.qs[i+j*m.ncatcg+1*scat*m.ncatcg][node.ID] = m.q1s[i+j*m.ncatcg]
				if node.Class == 1 {
					m.qs[i+j*m.ncatcg+2*scat*m.ncatcg][node.ID] = m.q2s[i+j*m.ncatcg]
					m.qs[i+j*m.ncatcg+3*scat*m.ncatcg][node.ID] = m.q2s[i+j*m.ncatcg]
				} else {
					m.qs[i+j*m.ncatcg+2*scat*m.ncatcg][node.ID] = m.q0s[i+j*m.ncatcg]
					m.qs[i+j*m.ncatcg+3*scat*m.ncatcg][node.ID] = m.q1s[i+j*m.ncatcg]
				}
			}
		}
	}
}

// updateProportions updates proportions if model parameters are
// changing.
func (m *BranchSiteGamma) updateProportions() {
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	bothcat := m.ncatcg * scat
	p0 := m.p0prop * m.p01sum
	p1 := m.p01sum - p0
	for c1 := 0; c1 < m.ncatsg; c1++ {
		for c2 := 0; c2 < m.ncatsg; c2++ {
			for c3 := 0; c3 < m.ncatsg; c3++ {
				for ecl := range m.gammac {
					catid := (((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl
					pq := m.gammasprop[c1] * m.gammasprop[c2] * m.gammasprop[c3] * m.gammacprop[ecl]
					m.prop[0][catid+bothcat*0] = pq * p0
					m.prop[0][catid+bothcat*1] = pq * p1
					m.prop[0][catid+bothcat*2] = pq * (1 - p0 - p1) * p0 / (p0 + p1)
					m.prop[0][catid+bothcat*3] = pq * (1 - p0 - p1) * p1 / (p0 + p1)
				}
			}
		}
	}

	for _, node := range m.data.Tree.NodeIDArray() {
		if node == nil {
			continue
		}
		scale := 0.0
		for i := 0; i < bothcat*4; i++ {
			scale += m.prop[0][i] * m.qs[i][node.ID].Scale
		}
		m.scale[node.ID] = scale
	}
	m.propdone = true
	m.expAllBr = false
}

// fillMatrices sets Q-matrices for all the rates.
func (m *BranchSiteGamma) fillMatrices(omega float64, dest []*codon.EMatrix) {
	for c1 := 0; c1 < m.ncatsg; c1++ {
		m.tmp[0] = m.gammas[c1]
		for c2 := 0; c2 < m.ncatsg; c2++ {
			m.tmp[1] = m.gammas[c2]
			for c3 := 0; c3 < m.ncatsg; c3++ {
				m.tmp[2] = m.gammas[c3]

				e := codon.NewEMatrix(m.data.cFreq)
				Q, s := codon.CreateRateTransitionMatrix(m.data.cFreq, m.kappa, omega, m.tmp, e.Q)
				e.Set(Q, s)
				err := e.Eigen()
				if err != nil {
					log.Fatal(err)
				}

				for ecl, rate := range m.gammac {
					catid := (((c1*m.ncatsg)+c2)*m.ncatsg+c3)*m.ncatcg + ecl

					e.Copy(dest[catid])
					dest[catid].ScaleD(rate)
				}
			}
		}

	}
}

// updateMatrices updates matrices if model parameters are changing.
func (m *BranchSiteGamma) updateMatrices() {
	if !m.q0done {
		m.fillMatrices(m.omega0, m.q0s)
		m.q0done = true
	}

	if !m.q1done {
		m.fillMatrices(1, m.q1s)
		m.q1done = true
	}

	if !m.q2done {
		m.fillMatrices(m.omega2, m.q2s)
		m.q2done = true
	}

	m.propdone = false
	m.expAllBr = false
}

// computePropBEB computes class proportions given the i, j and d parameters.
func (m *BranchSiteGamma) computePropBEB(prop []float64, i, j, d int) []float64 {
	if prop == nil {
		nClass := m.GetNClass()
		prop = make([]float64, nClass)
	}
	p0 := (1. + float64(j)/2*3 + float64(j%2)) / (3 * float64(d))
	p1 := (1. + float64(d-1-i)*3 + float64(j%2)) / (3 * float64(d))
	p2a := (1 - p0 - p1) * p0 / (p0 + p1)
	p2b := (1 - p0 - p1) * p1 / (p0 + p1)
	scat := m.ncatsg * m.ncatsg * m.ncatsg
	bothcat := m.ncatcg * scat
	for i := 0; i < bothcat; i++ {
		prop[i+bothcat*0] = p0 / float64(bothcat)
		prop[i+bothcat*1] = p1 / float64(bothcat)
		prop[i+bothcat*2] = p2a / float64(bothcat)
		prop[i+bothcat*3] = p2b / float64(bothcat)
	}

	return prop
}

//siteLMatrix returns site likelihood matrix. First index is w0, second
//w2, third class, forth position.
func (m *BranchSiteGamma) siteLMatrix(w0, w2 []float64) (res [][][][]float64) {
	res = make([][][][]float64, len(w0))
	nClass := m.GetNClass()
	nPos := m.data.cSeqs.Length()
	nni := m.data.Tree.MaxNodeID() + 1

	scat := m.ncatsg * m.ncatsg * m.ncatsg
	bothcat := m.ncatcg * scat

	// temporary storage for likelihood computation
	plh := make([][]float64, nni)
	for i := 0; i < nni; i++ {
		plh[i] = make([]float64, m.data.cFreq.GCode.NCodon+1)
	}

	counter := 0

	nWorkers := runtime.GOMAXPROCS(0)
	done := make(chan struct{}, nWorkers)
	type bebtask struct {
		iW0, iW2, class, pos int
	}
	tasks := make(chan bebtask, nPos)

	for i := 0; i < nWorkers; i++ {
		go func() {
			nni := m.data.Tree.MaxNodeID() + 1
			plh := make([][]float64, nni)
			for i := 0; i < nni; i++ {
				plh[i] = make([]float64, m.data.cFreq.GCode.NCodon+1)
			}
			for task := range tasks {
				res[task.iW0][task.iW2][task.class][task.pos] = m.fullSubL(task.class, task.pos, plh)
				done <- struct{}{}
			}
		}()
	}

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
				case class/bothcat == 0 && iW2 != 0:
					res[iW0][iW2][class] = res[iW0][0][class]
				case class/bothcat == 1 && (iW0 != 0 || iW2 != 0):
					res[iW0][iW2][class] = res[0][0][class]
				case class/bothcat == 3 && iW0 != 0:
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

// BEBPosterior returns BEB posterior values.
func (m *BranchSiteGamma) BEBPosterior() (res []float64) {
	nPos := m.data.cSeqs.Length()
	nClass := m.GetNClass()
	res = make([]float64, nPos)

	scat := m.ncatsg * m.ncatsg * m.ncatsg
	bothcat := m.ncatcg * scat

	log.Info("w0 and w2 grid")

	// first create a grid of parameters
	w0s := floatRange(0.05, 0.1, 10)
	log.Infof("w0: %v", strFltSlice(w0s))

	w2s := floatRange(1.5, 1, 10)
	log.Infof("w2: %v", strFltSlice(w2s))

	// we need to compute all the scaling
	// factors first
	m.expBranchesIfNeeded()
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

						s := 0.0
						for i := 0; i < bothcat; i++ {
							cl := i + bothcat*2
							s += matr[iW0][iW2][cl][pos] * prop[cl]
							cl = i + bothcat*3
							s += matr[iW0][iW2][cl][pos] * prop[cl]
						}
						tmp.SetFloat64(s)
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

// cGammaRates returns array with codon gamma rates using NEB approach.
func (m *BranchSiteGamma) cGammaRates() []float64 {
	clRates := make([]float64, m.GetNClass())
	scat := m.ncatsg * m.ncatsg * m.ncatsg

	for i := 0; i < m.ncatcg; i++ {
		for j := 0; j < scat; j++ {
			for k := 0; k < 4; k++ {
				clRates[m.ncatcg*scat*k+j*m.ncatcg+i] = m.gammac[i]
			}
		}
	}

	return m.NEBPosterior(clRates)
}

// Final prints NEB results (only if with positive selection).
func (m *BranchSiteGamma) Final(neb, beb, codonRates, siteRates, codonOmega bool) {
	startTime := time.Now()
	defer func() { m.summary.PosteriorTime = time.Since(startTime).Seconds() }()

	if m.ncatcg > 1 && codonRates {
		m.summary.CodonGammaRates = m.cGammaRates()
		log.Notice("Codon gamma rates posterior:", m.summary.CodonGammaRates)
	}

	if (neb || beb) && !m.fixw2 {
		classes := make([]float64, m.GetNClass())
		scat := m.ncatsg * m.ncatsg * m.ncatsg
		for i := 0; i < m.ncatcg; i++ {
			for j := 0; j < scat; j++ {
				classes[m.ncatcg*scat*2+i+j*m.ncatcg] = 1
				classes[m.ncatcg*scat*3+i+j*m.ncatcg] = 1
			}
		}

		if neb {
			m.summary.SitePosteriorNEB = m.NEBPosterior(classes)

			log.Notice("NEB analysis")
			m.PrintPosterior(m.summary.SitePosteriorNEB)
		}

		if beb {
			m.summary.SitePosteriorBEB = m.BEBPosterior()

			log.Notice("BEB analysis")
			m.PrintPosterior(m.summary.SitePosteriorBEB)
		}
	}
	if m.ncatcg > 1 && m.proportional {
		m.update()
		log.Infof("Codon rates: %v", m.gammac)
		log.Infof("Codon props: %v", m.gammacprop)
	}
	if m.ncatsg > 1 && m.proportional {
		m.update()
		log.Infof("Site rates: %v", m.gammas)
		log.Infof("Site props: %v", m.gammasprop)
	}
}

// update updates matrices and proportions.
func (m *BranchSiteGamma) update() {
	if !m.gammasdone {
		if m.ncatsg > 1 {
			if !m.proportional {
				m.gammas = dist.DiscreteGamma(m.alphas, m.alphas, m.ncatsg, false, m.tmp, m.gammas)
			} else {
				rs1s := math.Exp(m.rs1s)
				rs2s := math.Exp(m.rs2s)

				synScale := rs1s*m.ps1s + 1*(1-m.ps1s)*m.ps2s + rs2s*(1-m.ps1s)*(1-m.ps2s)
				m.gammas[0] = rs1s / synScale
				m.gammas[1] = 1 / synScale
				m.gammas[2] = rs2s / synScale
				m.gammasprop[0] = m.ps1s
				m.gammasprop[1] = (1 - m.ps1s) * m.ps2s
				m.gammasprop[2] = (1 - m.ps1s) * (1 - m.ps2s)
			}
			m.q0done = false
			m.q1done = false
			m.q2done = false
		} else {
			m.gammas[0] = 1
			m.gammasprop[0] = 1
		}
		m.gammasdone = true
	}
	if !m.gammacdone {
		if m.ncatcg > 1 {
			if !m.proportional {
				m.gammac = dist.DiscreteGamma(m.alphac, m.alphac, m.ncatcg, false, m.tmp, m.gammac)
			} else {
				rs1c := math.Exp(m.rs1c)
				rs2c := math.Exp(m.rs2c)

				synScale := rs1c*m.ps1c + 1*(1-m.ps1c)*m.ps2c + rs2c*(1-m.ps1c)*(1-m.ps2c)
				m.gammac[0] = rs1c / synScale
				m.gammac[1] = 1 / synScale
				m.gammac[2] = rs2c / synScale
				m.gammacprop[0] = m.ps1c
				m.gammacprop[1] = (1 - m.ps1c) * m.ps2c
				m.gammacprop[2] = (1 - m.ps1c) * (1 - m.ps2c)
			}
			m.q0done = false
			m.q1done = false
			m.q2done = false
		} else {
			m.gammac[0] = 1
			m.gammacprop[0] = 1
		}
		m.gammacdone = true
	}
	if !m.q0done || !m.q1done || !m.q2done {
		m.updateMatrices()
	}
	if !m.propdone {
		m.updateProportions()
	}
}

// Summary returns the run summary (site posterior for NEB and BEB).
func (m *BranchSiteGamma) Summary() interface{} {
	if !m.summary.Empty() {
		return m.summary
	}
	// nil prevents json from printing "{}"
	return nil
}

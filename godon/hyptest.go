package main

import (
	"bitbucket.org/Davydov/godon/cmodel"
)

const (
	// defaultSThr  is the minimal D-statistics value, for which
	// beb should be computed (qchisq(0.9, df=1)).
	defaultSThr = "2.705543"
)

// hypTest performs the hypthesis testing
func hypTest() (tests []HypTestSummary, optimizations []OptimizationSummary) {
	//transfer options from hypTest command
	alignmentFileName = hTestAlignmentFileName
	treeFileName = hTestTreeFileName
	model = hTestModel

	data, err := newData()
	if err != nil {
		log.Fatal(err)
	}

	if *m0Tree && !*noOptBrLen {
		m0ms := newModelSettings(data)
		if *m0TreeGamma {
			m0ms.name = "M0G"
		} else {
			m0ms.name = "M0"
		}
		m0ms.startF = ""
		m0model, err := m0ms.createInitalized(false)
		if err != nil {
			log.Fatal(err)
		}
		m0opt := newOptimizerSettings(m0model)
		log.Notice("Optimizing branch lengths using M0")
		res := runOptimization(m0model, m0opt, nil, true)
		optimizations = append(optimizations, res)
		*noOptBrLen = true
	}

	if (*model == "BS" || *model == "BSG") && (data.GetNClass1() < 1 || *testAllBranches) {
		log.Notice("Testing multiple branches")
		toTest := make([]int, 0, data.Tree.NNodes())

		// first compute branches to test; and set all branches to class 0
		for node := range data.Tree.Walker(nil) {
			node.Class = 0
			if node.IsRoot() {
				continue
			}
			if *noLeavesTest && node.IsTerminal() {
				continue
			}
			toTest = append(toTest, node.ID)
		}

		nodes := data.Tree.NodeIDArray()

		for i, nid := range toTest {
			log.Noticef("Testing branch %d/%d", i+1, len(toTest))
			nodes[nid].Class = 1
			log.Noticef("Foreground branch: %s", data.Tree.ShortClassString())
			tests = append(tests, performSingleTest(data))
			nodes[nid].Class = 0
		}
	} else {
		tests = append(tests, performSingleTest(data))
	}
	return tests, optimizations
}

// performSingleTest preforms a test for given data
func performSingleTest(data *cmodel.Data) (summary HypTestSummary) {
	summary.Tree = data.Tree.ClassString()

	// name of the H1 extra parameter
	var extraPar string
	switch {
	case *model == "BS" || *model == "BSG":
		extraPar = "omega2"
	case *model == "M8":
		extraPar = "omega"
	default:
		log.Fatalf("Unknown model '%v'", *model)
	}

	ms := newModelSettings(data)

	ms.fixw = true
	m0, err := ms.createInitalized(true)
	if err != nil {
		log.Fatal(err)
	}

	o0 := newOptimizerSettings(m0)

	log.Notice("Running H0")
	res0 := runOptimization(m0, o0, nil, true)
	res0.Hypothesis = "H0"
	summary.Optimizations = append(summary.Optimizations, res0)

	ms.fixw = false
	m1, err := ms.createInitalized(true)
	if err != nil {
		log.Fatal(err)
	}
	o1 := newOptimizerSettings(m1)

	log.Notice("Running H1")
	res1 := runOptimization(m1, o1, nil, true)
	res1.Hypothesis = "H1"
	summary.Optimizations = append(summary.Optimizations, res1)

	var l0, l1 float64
	l0 = res0.Optimizer.GetMaxLikelihood()
	l1 = res1.Optimizer.GetMaxLikelihood()

	log.Noticef("Starting with D=%g", 2*(l1-l0))

	if *quick {
		o0.method = "none"
		o1.method = "none"
	}

	for updated, justUpdatedH0 := true, false; updated; {
		// stop if nothing has been updated.
		// this loop normally should run for a single time.
		updated = false

		// if l1 < l0, rerun H1 starting from H0.
		if lrt := 2 * (l1 - l0); lrt < 0 {
			updated = true
			justUpdatedH0 = false

			h0par := res0.Optimizer.GetMaxLikelihoodParameters()
			h0par[extraPar] = 1

			log.Noticef("Rerunning H1 because of negative LR (D=%g)",
				lrt)
			res1 = runOptimization(m1, o1, h0par, true)
			res1.Hypothesis = "H1"
			summary.Optimizations = append(summary.Optimizations, res1)
			l1 = res1.Optimizer.GetMaxLikelihood()
			// some optimizers (e.g. n_lbfgsb) can return a maxL
			// which is worse than starting point
			if newLrt := 2 * (l1 - l0); newLrt < 0  {
				if o1.method != "none" {
					log.Warning("Warning: optimizer failed to correct negative LR, fallback to the quick mode")
					o1.method = "none"
				} else {
					log.Warning("Warning: couldn't get rid of negative LR; models are not strictly nested")
					break
				}
			}
		}

		// if significant (D>thr), rerun H0 starting from H1
		if lrt := 2 * (l1 - l0); lrt > *sThr && !*quick && !justUpdatedH0 {
			// prevent multiple updates of H0 starting from the same
			// H1-like point
			justUpdatedH0 = true

			h1par := res1.Optimizer.GetMaxLikelihoodParameters()
			delete(h1par, extraPar)
			log.Noticef("Rerunning H0, trying to reduce LR (D=%g)",
				lrt)
			res0Alt := runOptimization(m0, o0, h1par, true)
			res0Alt.Hypothesis = "H0"
			summary.Optimizations = append(summary.Optimizations, res0Alt)
			l0Alt := res0Alt.Optimizer.GetMaxLikelihood()
			if l0Alt > l0 {
				updated = true
				res0 = res0Alt
				l0 = l0Alt
			}
		}
	}

	// get rid of sligtly positive LRT; this should require maximum one extra
	// likelihood computation
	if lrt := 2 * (l1 - l0); lrt > 0 && lrt <= *sThr {
		h1par := res1.Optimizer.GetMaxLikelihoodParameters()
		delete(h1par, extraPar)
		o0.method = "none"
		log.Noticef("Rerunning H0, trying to reduce small LR (D=%g)",
			lrt)
		res0Alt := runOptimization(m0, o0, h1par, true)
		res0Alt.Hypothesis = "H0"
		summary.Optimizations = append(summary.Optimizations, res0Alt)
		l0Alt := res0Alt.Optimizer.GetMaxLikelihood()
		if l0Alt > l0 {
			res0 = res0Alt
			l0 = l0Alt
		}
	}

	// one last round of getting rid of negative lrt, maximum one extra
	// likelihood computation
	if lrt := 2 * (l1 - l0); lrt < 0 {
		h0par := res0.Optimizer.GetMaxLikelihoodParameters()
		h0par[extraPar] = 1
		o1.method = "none"

		log.Noticef("Rerunning H1 because of negative LR (D=%g)",
			lrt)
		res1 = runOptimization(m1, o1, h0par, true)
		res1.Hypothesis = "H1"
		summary.Optimizations = append(summary.Optimizations, res1)
		l1 = res1.Optimizer.GetMaxLikelihood()
	}

	log.Noticef("Final D=%g", 2*(l1-l0))

	// final stores BEB & NEB results
	var final0Summary, final1Summary interface{}

	if 2*(l1-l0) <= *sThr {
		// not significant, no need for NEB & BEB
		*neb = false
		*beb = false
	}

	// final BEB & NEB other computations
	if *final {
		h0par := res0.Optimizer.GetMaxLikelihoodParameters()
		final0Summary = computeFinal(m0, h0par)

		h1par := res1.Optimizer.GetMaxLikelihoodParameters()
		final1Summary = computeFinal(m1, h1par)
	}

	summary.H0 = HypSummary{
		MaxLnL:         res0.Optimizer.GetMaxLikelihood(),
		MaxLParameters: res0.Optimizer.GetMaxLikelihoodParameters(),
		Final:          final0Summary,
	}
	summary.H1 = HypSummary{
		MaxLnL:         res1.Optimizer.GetMaxLikelihood(),
		MaxLParameters: res1.Optimizer.GetMaxLikelihoodParameters(),
		Final:          final1Summary,
	}

	log.Noticef("lnL0=%f, lnL1=%f",
		summary.H0.MaxLnL,
		summary.H1.MaxLnL)

	return
}

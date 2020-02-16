package main

import (
	"encoding/json"

	"bitbucket.org/Davydov/godon/checkpoint"
	"bitbucket.org/Davydov/godon/cmodel"
)

const (
	// defaultSThr  is the minimal D-statistics value, for which
	// beb should be computed (qchisq(0.9, df=1)).
	defaultSThr = "2.705543"

	// minLrt is the minimal tolerated LRT value.
	minLrt = 1e-6

	// maxOptimizations is the maximum number of likelihood
	// optimizations to perform
	maxOptimizations = 16
)

// thoroughMin sets minimal starting value likelihood considering
// thorough flag.
func thoroughMin(min float64) float64 {
	if *thorough {
		return 1
	}
	return min
}

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

	if (*m0Tree || *m0TreeGamma) && !*noOptBrLen {
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
		key := m0ms.name + ":" + data.Tree.ShortClassString()
		res := runOptimization(m0model, m0opt, nil, 1, []byte(key), true)
		optimizations = append(optimizations, res)
		*noOptBrLen = true
	}

	if (*model == "BS" || *model == "BSG") && (data.GetNClass1() < 1 || *testAllBranches) {
		// set testAllBranches so we know multiple branches were tested
		*testAllBranches = true
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

		if len(toTest) == 0 {
			log.Warningf("No branches to test")
		}
	} else {
		// update testAllBranches so we know that only one test was performed
		*testAllBranches = false
		tests = append(tests, performSingleTest(data))
	}
	return tests, optimizations
}

// saveSummary saves summary to checkpoint
func saveSummary(summary interface{}, key []byte) {
	var b []byte = []byte{}
	var err error
	if summary != nil {
		b, err = json.Marshal(summary)
	}
	if err != nil {
		log.Error("Error marshalling final summary", err)
	} else {
		err = checkpoint.SaveData(checkpointDB, key, b)
		if err != nil {
			log.Error("Error saving final summary", err)
		}
	}
}

// loadSummary loads summary pseudoobjects
func loadSummary(keyFinal0, keyFinal1 []byte) (*checkpoint.PseudoObject, *checkpoint.PseudoObject, bool) {
	final0SummaryB, _ := checkpoint.LoadData(checkpointDB, keyFinal0)
	final1SummaryB, _ := checkpoint.LoadData(checkpointDB, keyFinal1)
	if final0SummaryB == nil || final1SummaryB == nil {
		return nil, nil, false
	}
	final0Summary := checkpoint.NewPseudoObject(final0SummaryB)
	final1Summary := checkpoint.NewPseudoObject(final1SummaryB)
	return final0Summary, final1Summary, true
}

// performSingleTest preforms a test for given data
func performSingleTest(data *cmodel.Data) (summary HypTestSummary) {
	summary.Tree = data.Tree.ClassString()
	clstr := data.Tree.ShortClassString()

	keyFinal0 := []byte(*model + ":" + "H0" + ":final:" + clstr)
	keyFinal1 := []byte(*model + ":" + "H1" + ":final:" + clstr)

	final0SummarySaved, final1SummarySaved, finalAvail := loadSummary(keyFinal0, keyFinal1)

	// names and default values of the H1 extra parameters
	extraPar := make(map[string]float64, 2)
	switch {
	case *model == "BS" || *model == "BSG":
		extraPar["omega2"] = 1
	case *model == "M8":
		extraPar["omega"] = 1
	case *model == "M2a":
		extraPar["omega2"] = 1
		extraPar["p1prop"] = 1
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
	if finalAvail {
		o0.method = "none"
	}

	log.Notice("Running H0")
	key0 := *model + ":" + "H0" + ":" + clstr
	res0 := runOptimization(m0, o0, nil, 1, []byte(key0), true)
	res0.Hypothesis = "H0"
	summary.Optimizations = append(summary.Optimizations, res0)

	ms.fixw = false
	m1, err := ms.createInitalized(true)
	if err != nil {
		log.Fatal(err)
	}
	o1 := newOptimizerSettings(m1)
	if finalAvail {
		o1.method = "none"
	}

	log.Notice("Running H1")
	key1 := *model + ":" + "H1" + ":" + clstr
	res1 := runOptimization(m1, o1, nil, 1, []byte(key1), true)
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

	for updated, justUpdatedH0 := true, false; updated && !finalAvail; {
		// stop if nothing has been updated.
		// this loop normally should run for a single time.
		updated = false

		// if l1 < l0, rerun H1 starting from H0.
		if lrt := 2 * (l1 - l0); lrt < 0 {
			updated = true
			justUpdatedH0 = false

			h0par := res0.Optimizer.GetMaxLikelihoodParameters()
			for parName, parVal := range extraPar {
				h0par[parName] = parVal
			}

			log.Noticef("Rerunning H1 because of negative LR (D=%g)",
				lrt)
			res1 = runOptimization(m1, o1, h0par, thoroughMin(l1), nil, true)
			res1.Hypothesis = "H1"
			summary.Optimizations = append(summary.Optimizations, res1)
			l1 = res1.Optimizer.GetMaxLikelihood()
			// some optimizers (e.g. n_lbfgsb) can return a maxL
			// which is worse than starting point
			if newLrt := 2 * (l1 - l0); newLrt < 0 {
				if o1.method != "none" {
					log.Warning("Warning: optimizer failed to correct negative LR, fallback to the quick mode")
					o1.method = "none"
				} else {
					log.Warning("Warning: couldn't get rid of negative LR; models are not strictly nested")
					break
				}
			}
		}

		// stop if too many updates; we prefer to have slightly positive LRT
		if len(summary.Optimizations) >= maxOptimizations {
			log.Warningf("Performed %s optimization, finishing")
			break
		}

		// if significant (D>thr), rerun H0 starting from H1
		if lrt := 2 * (l1 - l0); (lrt > *sThr || (lrt > minLrt && !*thorough)) && !justUpdatedH0 {
			// prevent multiple updates of H0 starting from the same
			// H1-like point
			justUpdatedH0 = true

			h1par := res1.Optimizer.GetMaxLikelihoodParameters()
			for parName := range extraPar {
				delete(h1par, parName)
			}

			log.Noticef("Rerunning H0, trying to reduce LR (D=%g)",
				lrt)
			res0Alt := runOptimization(m0, o0, h1par, thoroughMin(l0), nil, true)
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
	if lrt := 2 * (l1 - l0); lrt > minLrt && (lrt <= *sThr || !*thorough) && !finalAvail {
		h1par := res1.Optimizer.GetMaxLikelihoodParameters()
		for parName := range extraPar {
			delete(h1par, parName)
		}

		o0.method = "none"
		log.Noticef("Rerunning H0, trying to reduce LR (D=%g)",
			lrt)
		res0Alt := runOptimization(m0, o0, h1par, thoroughMin(l0), nil, true)
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
	if lrt := 2 * (l1 - l0); lrt < 0 && !finalAvail {
		h0par := res0.Optimizer.GetMaxLikelihoodParameters()
		for parName, parVal := range extraPar {
			h0par[parName] = parVal
		}

		o1.method = "none"

		log.Noticef("Rerunning H1 because of negative LR (D=%g)",
			lrt)
		res1 = runOptimization(m1, o1, h0par, thoroughMin(l1), nil, true)
		res1.Hypothesis = "H1"
		summary.Optimizations = append(summary.Optimizations, res1)
		l1 = res1.Optimizer.GetMaxLikelihood()
	}

	// checkpoint parameters after all the optimizations
	h0par := res0.Optimizer.GetMaxLikelihoodParameters()
	checkpointParameters(m0, o0, h0par, []byte(key0))
	h1par := res1.Optimizer.GetMaxLikelihoodParameters()
	checkpointParameters(m1, o1, h1par, []byte(key1))

	log.Noticef("Final D=%g", 2*(l1-l0))

	// final stores BEB & NEB results
	var final0Summary, final1Summary interface{}

	// should we compute BEB & NEB for this run?
	runNEB := *neb
	runBEB := *beb

	if 2*(l1-l0) <= *sThr {
		// not significant, no need for NEB & BEB
		runNEB = false
		runBEB = false
	}

	// final BEB & NEB other computations
	if *final {
		if finalAvail {
			final0Summary = final0SummarySaved.SelfOrNil()
			final1Summary = final1SummarySaved.SelfOrNil()
		} else {
			final0Summary = computeFinal(m0, h0par, runNEB, runBEB)
			saveSummary(final0Summary, keyFinal0)
			final1Summary = computeFinal(m1, h1par, runNEB, runBEB)
			saveSummary(final1Summary, keyFinal1)
		}
	} else {
		// even if no finalSummary is needed, we need to indicate that
		// all optimizations were finished
		saveSummary(nil, keyFinal0)
		saveSummary(nil, keyFinal1)
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

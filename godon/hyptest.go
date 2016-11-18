package main

const (
	// defaultSThr  is the minimal D-statistics value, for which
	// beb should be computed (qchisq(0.9, df=1)).
	defaultSThr = "2.705543"
)

// hypTest performs the hypthesis testing
func hypTest() (summary HypTestSummary) {
	//transfer options from hypTest command
	alignmentFileName = hTestAlignmentFileName
	treeFileName = hTestTreeFileName
	model = hTestModel

	// don't compute BEB & NEB
	globalNoFinal := *noFinal
	*noFinal = true

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

	log.Notice("Running H0")
	res0 := runOptimization(true, nil)
	res0.Hypothesis = "H0"
	summary.Runs = append(summary.Runs, res0)

	log.Notice("Running H1")
	res1 := runOptimization(false, nil)
	res1.Hypothesis = "H1"
	summary.Runs = append(summary.Runs, res1)

	var l0, l1 float64
	l0 = res0.Optimizer.GetMaxLikelihood()
	l1 = res1.Optimizer.GetMaxLikelihood()

	for updated, justUpdatedH0 := true, false; updated; {
		// stop if nothing has been updated.
		// this loop normally should run for a single time.
		updated = false

		// if l1 < l0, rerun H1 starting from H0.
		if l1 < l0 {
			updated = true
			justUpdatedH0 = false

			h0par := res0.Optimizer.GetMaxLikelihoodParameters()
			h0par[extraPar] = 1

			if *quick {
				*method = "none"
			}
			log.Noticef("Rerunning H1 because of negative LRT (lnL1-lnL0=%f)",
				l1-l0)
			res1 = runOptimization(false, h0par)
			res1.Hypothesis = "H1"
			summary.Runs = append(summary.Runs, res1)
			l1 = res1.Optimizer.GetMaxLikelihood()
		}

		// if significant (D>thr), rerun H0 starting from H1
		if 2*(l1-l0) > *sThr && !*quick && !justUpdatedH0 {
			// prevent multiple updates of H0 starting from the same
			// H1-like point
			justUpdatedH0 = true

			h1par := res1.Optimizer.GetMaxLikelihoodParameters()
			delete(h1par, extraPar)
			log.Noticef("Rerunning H0, trying to reduce LRT (lnL1-lnL0=%f)",
				l1-l0)
			res0Alt := runOptimization(true, h1par)
			res0Alt.Hypothesis = "H0"
			summary.Runs = append(summary.Runs, res0Alt)
			l0Alt := res0Alt.Optimizer.GetMaxLikelihood()
			if l0Alt > l0 {
				updated = true
				res0 = res0Alt
				l0 = l0Alt
			}
		}
	}

	// final BEB & NEB computation
	if !globalNoFinal && 2*(l1-l0) > *sThr {
		h1par := res1.Optimizer.GetMaxLikelihoodParameters()
		*noFinal = false
		*method = "none"
		log.Notice("Computing posterior for H1")
		res1 = runOptimization(false, h1par)
		res1.Hypothesis = "H1"
		summary.Runs = append(summary.Runs, res1)
	}

	summary.H0 = HypSummary{
		MaxLnL:         res0.Optimizer.GetMaxLikelihood(),
		MaxLParameters: res0.Optimizer.GetMaxLikelihoodParameters(),
		Final:          res0.Model,
	}
	summary.H1 = HypSummary{
		MaxLnL:         res1.Optimizer.GetMaxLikelihood(),
		MaxLParameters: res1.Optimizer.GetMaxLikelihoodParameters(),
		Final:          res1.Model,
	}

	log.Noticef("lnL0=%f, lnL1=%f",
		summary.H0.MaxLnL,
		summary.H1.MaxLnL)

	return
}

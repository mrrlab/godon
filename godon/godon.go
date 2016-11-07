/*

Godon implements number of codon models (including M0 and
branch-site). It includes various likelihood optimizers as well as
Metropolis-Hastings sampler.

The basic usage of godon looks like this:

	godon alignment.fst tree.nwk

, this will run M0 model with a default optimizer (LBFGS-B).

You can change a model and an optimizer:

	godon -model BS -method simplex alignment.fst tree.nwk

The above will use the branch-site model and downhill simplex optimizer.

To see all the options run:

	godon -h

*/
package main

import (
	"bufio"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
	"time"

	"github.com/op/go-logging"

	"bitbucket.org/Davydov/godon/bio"
	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/tree"
)

// These three variables are set during the compilation.
var githash = ""
var gitbranch = ""
var buildstamp = ""
var version = fmt.Sprintf("branch: %s, revision: %s, build time: %s", gitbranch, githash, buildstamp)

// Logger settings.
var log = logging.MustGetLogger("godon")
var formatter = logging.MustStringFormatter(`%{message}`)

// lastLine returns the last line of a file content.
func lastLine(fn string) (line string, err error) {
	f, err := os.Open(fn)
	if err != nil {
		return line, err
	}
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line = scanner.Text()
	}
	err = scanner.Err()
	return line, err
}

// getAggModeFromString returns an aggregation mode constant from cmodel
// from a string.
func getAggModeFromString(aggModeString string) (cmodel.AggMode, error) {
	switch aggModeString {
	case "none":
		return cmodel.AggNone, nil
	case "observed":
		return cmodel.AggObserved, nil
	case "observed_new":
		return cmodel.AggObservedNew, nil
	case "fixed":
		return cmodel.AggFixed, nil
	case "random":
		return cmodel.AggRandom, nil
	}
	return cmodel.AggNone, fmt.Errorf("Unknown aggregation mode: %s", aggModeString)
}

// getModelFromString returns a model from string and other parameters.
func getModelFromString(model string, cali codon.Sequences, t *tree.Tree, cf codon.Frequency,
	fixw bool, ncatb, ncatsg, ncatcg int) (cmodel.TreeOptimizableSiteClass, error) {
	switch model {
	case "M0":
		log.Info("Using M0 model")
		return cmodel.NewM0(cali, t, cf), nil
	case "M1a":
		log.Info("Using M1a model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewM2(cali, t, cf, false, ncatsg, ncatcg), nil
	case "M2a":
		log.Info("Using M2a model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewM2(cali, t, cf, true, ncatsg, ncatcg), nil
	case "M7":
		log.Info("Using M7 model")
		log.Infof("%d beta categories, %d site gamma categories, %d codon gama categories", ncatb, ncatsg, ncatcg)
		return cmodel.NewM8(cali, t, cf, false, false, ncatb, ncatsg, ncatcg), nil
	case "M8":
		log.Info("Using M8 model")
		log.Infof("%d beta categories, %d site gamma categories, %d codon gama categories", ncatb, ncatsg, ncatcg)
		return cmodel.NewM8(cali, t, cf, true, fixw, ncatb, ncatsg, ncatcg), nil
	case "BSG":
		log.Info("Using branch site gamma model")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewBranchSiteGamma(cali, t, cf, fixw, ncatsg, ncatcg), nil
	case "BSGE":
		log.Info("Using branch site gamma model with explicit rates")
		log.Infof("%d site gamma categories, %d codon gama categories", ncatsg, ncatcg)
		return cmodel.NewBranchSiteGammaERates(cali, t, cf, fixw, ncatsg, ncatcg), nil
	case "BS":
		log.Info("Using branch site model")
		return cmodel.NewBranchSite(cali, t, cf, fixw), nil
	}
	return nil, errors.New("Unknown model specification")
}

// getOptimizerFromString returns an optimizer from a string.
func getOptimizerFromString(method string, accept, annealingSkip int, seed int64) (optimize.Optimizer, error) {
	switch method {
	case "lbfgsb":
		return optimize.NewLBFGSB(), nil
	case "simplex":
		return optimize.NewDS(), nil
	case "mh":
		chain := optimize.NewMH(false, 0)
		chain.AccPeriod = accept
		return chain, nil
	case "annealing":
		chain := optimize.NewMH(true, annealingSkip)
		chain.AccPeriod = accept
		return chain, nil
	case "n_lbfgs":
		return optimize.NewNLOPT(optimize.NLOPT_LBFGS, seed), nil
	case "n_simplex":
		return optimize.NewNLOPT(optimize.NLOPT_SIMPLEX, seed), nil
	case "n_cobyla":
		return optimize.NewNLOPT(optimize.NLOPT_COBYLA, seed), nil
	case "n_bobyqa":
		return optimize.NewNLOPT(optimize.NLOPT_BOBYQA, seed), nil
	case "n_sqp":
		return optimize.NewNLOPT(optimize.NLOPT_SQP, seed), nil
	case "n_direct":
		return optimize.NewNLOPT(optimize.NLOPT_DIRECT, seed), nil
	case "n_crs":
		return optimize.NewNLOPT(optimize.NLOPT_CRS, seed), nil
	case "n_mlsl":
		return optimize.NewNLOPT(optimize.NLOPT_MLSL, seed), nil
	case "none":
		return optimize.NewNone(), nil
	}
	return nil, fmt.Errorf("Unknown optimization method: %s", method)
}

func main() {
	startTime := time.Now()

	summary := Summary{}

	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Godon %s\n", version)
		fmt.Fprintf(os.Stderr, "Usage:\n")
		flag.PrintDefaults()
	}

	// model
	model := flag.String("model", "M0", "todel type (M0 or BS for branch site)")
	gcodeID := flag.Int("gcode", 1, "NCBI genetic code id, standard by default")
	fgBranch := flag.Int("fg", -1, "fg branch number")
	maxBrLen := flag.Float64("maxbrlen", 100, "maximum branch length")
	noOptBrLen := flag.Bool("nobrlen", false, "don't optimize branch lengths")
	cFreq := flag.String("cfreq", "F3X4", "codon frequecny (F0 or F3X4)")
	cFreqFileName := flag.String("cfreqfn", "", "codon frequencies file (overrides -cfreq)")
	fixw := flag.Bool("fixw", false, "fix omega=1 (for the branch-site and M8 models)")
	ncatsg := flag.Int("ncatsg", 1, "number of categories for the site gamma rate variation (no variation by default)")
	ncatcg := flag.Int("ncatcg", 1, "number of categories for the codon gamma rate variation (no variation by default)")
	ncatb := flag.Int("ncatb", 4, "number of the categories for the beta distribution (models M7&M8)")
	noFinal := flag.Bool("nofinal", false, "don't perform final extra computations, i.e. NEB and BEB site posterior")

	// optimizer parameters
	randomize := flag.Bool("randomize", false, "use uniformly distributed random starting point; "+
		"by default random starting point is distributed around realistic parameter values")
	iterations := flag.Int("iter", 10000, "number of iterations")
	report := flag.Int("report", 10, "report every N iterations")
	method := flag.String("method", "lbfgsb", "optimization method to use "+
		"(lbfgsb: limited-memory Broyden–Fletcher–Goldfarb–Shanno with bounding constraints, "+
		"simplex: downhill simplex, "+
		"annealing: simullated annealing, "+
		"mh: Metropolis-Hastings, "+
		"n_lbfgs: LBFGS from nlopt, "+
		"n_simplex: downhill simplex from nlopt, "+
		"n_cobyla: COBYLA from nlopt, "+
		"n_bobyqa: BOBYQA from nlopt, "+
		"n_sqp: SQP from nlopt, "+
		"n_mlsl: MLSL from nlopt (BOBYQA local optimizer), "+
		"none: just compute likelihood, no optimization"+
		"); you can chain optimizers with a plus sign")
	// mcmc parameters
	accept := flag.Int("accept", 200, "report acceptance rate every N iterations")

	// adaptive mcmc parameters
	adaptive := flag.Bool("adaptive", false, "use adaptive MCMC")
	skip := flag.Int("skip", -1, "number of iterations to skip for adaptive mcmc (5% by default)")
	maxAdapt := flag.Int("maxadapt", -1, "stop adapting after iteration (20% by default)")

	// optimizations
	aggregate := flag.String("aggregate", "none", "state aggregation mode: "+
		"observed (all positions, keep observed states), "+
		"observed_new (new implementation of observed), "+
		"fixed (absolutely conserved positions, keep observed), "+
		"random (like observed, but non-aggregated states are shuffled between the positions)")
	printFull := flag.Bool("printfull", false, "print full (non-aggregated) likelihood in the end of optimization")

	// technical
	nThreads := flag.Int("nt", 0, "number of threads to use")
	seed := flag.Int64("seed", -1, "random generator seed, default time based")
	cpuProfile := flag.String("cpuprofile", "", "write cpu profile to file")

	// input/output
	outLogF := flag.String("log", "", "write log to a file")
	outF := flag.String("out", "", "write optimization trajectory to a file")
	outTreeF := flag.String("tree", "", "write tree to a file")
	startF := flag.String("start", "", "read start position from the trajectory or JSON file")
	logLevel := flag.String("loglevel", "notice", "loglevel ('critical', 'error', 'warning', 'notice', 'info', 'debug')")
	jsonF := flag.String("json", "", "write json output to a file")

	flag.Parse()

	// logging
	logging.SetFormatter(formatter)

	var backend *logging.LogBackend
	if *outLogF != "" {
		f, err := os.OpenFile(*outLogF, os.O_WRONLY|os.O_CREATE|os.O_APPEND, 0666)
		if err != nil {
			log.Fatal("Error creating log file:", err)
		}
		defer f.Close()
		backend = logging.NewLogBackend(f, "", 0)
	} else {
		backend = logging.NewLogBackend(os.Stderr, "", 0)
	}
	logging.SetBackend(backend)

	level, err := logging.LogLevel(*logLevel)
	if err != nil {
		log.Fatal(err)
	}
	logging.SetLevel(level, "godon")
	logging.SetLevel(level, "optimize")
	logging.SetLevel(level, "cmodel")

	// print revision
	log.Info(version)
	summary.Version = version

	// print commandline
	log.Info("Command line:", os.Args)
	summary.CommandLine = os.Args

	if *seed == -1 {
		*seed = time.Now().UnixNano()
		log.Debug("Random seed from time")
	}
	log.Infof("Random seed=%v", *seed)
	summary.Seed = *seed

	rand.Seed(*seed)
	runtime.GOMAXPROCS(*nThreads)

	effectiveNThreads := runtime.GOMAXPROCS(0)
	log.Infof("Using threads: %d.\n", effectiveNThreads)
	summary.NThreads = effectiveNThreads

	if len(flag.Args()) < 2 {
		log.Fatal("you should specify tree and alignment")
		return
	}

	gcode, ok := bio.GeneticCodes[*gcodeID]
	if !ok {
		log.Fatalf("couldn't load genetic code with id=%d", gcodeID)
	}
	log.Infof("Genetic code: %d, \"%s\"", gcode.ID, gcode.Name)

	fastaFile, err := os.Open(flag.Args()[0])
	if err != nil {
		log.Fatal(err)
	}
	defer fastaFile.Close()

	ali, err := bio.ParseFasta(fastaFile)
	if err != nil {
		log.Fatal(err)
	}

	cali, err := codon.ToCodonSequences(ali, gcode)
	if err != nil {
		log.Fatal(err)
	}

	if cali.Length() == 0 {
		log.Fatal("Zero length alignment")
	}
	log.Infof("Read alignment of %d codons, %d fixed positions, %d ambiguous positions", cali.Length(), cali.NFixed(), cali.NAmbiguous())

	treeFile, err := os.Open(flag.Args()[1])
	if err != nil {
		log.Fatal(err)
	}
	defer treeFile.Close()

	t, err := tree.ParseNewick(treeFile)
	if err != nil {
		log.Fatal(err)
	}

	log.Debugf("intree=%s", t)
	log.Debugf("brtree=%s", t.BrString())

	// Root the tree in the end
	var root = false
	var rootID = 0

	if t.IsRooted() {
		log.Warning("Tree is rooted. Will unroot.")
		root = true
		rootID, err = t.Unroot()
		if err != nil {
			log.Fatal("Error unrooting tree:", err)
		}
	}

	log.Infof("intree_unroot=%s", t)
	summary.StartTree = t.String()
	log.Debugf("brtree_unroot=%s", t.BrString())
	log.Debug(t.FullString())

	var cf codon.Frequency

	if *cFreqFileName != "" {
		cFreqFile, err := os.Open(*cFreqFileName)
		if err != nil {
			log.Fatal(err)
		}
		cf, err = codon.ReadFrequency(cFreqFile, gcode)
		if err != nil {
			log.Fatal(err)
		}
	} else {
		switch *cFreq {
		case "F0":
			log.Info("F0 frequency")
			cf = codon.F0(cali)
		case "F3X4":
			log.Info("F3X4 frequency")
			cf = codon.F3X4(cali)
		default:
			log.Fatal("Unknow codon freuquency specification")
		}
	}
	log.Debug(cf)

	if *fgBranch >= 0 {
		for _, node := range t.NodeIDArray() {
			if node == nil {
				continue
			}
			if node.ID == *fgBranch {
				node.Class = 1
			} else {
				node.Class = 0
			}
		}
	} else {
		class1 := 0
		for range t.ClassNodes(1) {
			class1++
		}
		if class1 == 0 &&
			(*model == "BS" || *model == "BSG" || *model == "BSGE") {
			log.Warning("Warning: no class=1 nodes")
		}
	}

	if *cpuProfile != "" {
		f, err := os.Create(*cpuProfile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	m, err := getModelFromString(*model, cali, t, cf, *fixw, *ncatb, *ncatsg, *ncatcg)
	if err != nil {
		log.Fatal(err)
	}

	log.Infof("Model has %d site class(es)", m.GetNClass())

	if !*noOptBrLen {
		log.Info("Will optimize branch lengths")
		log.Infof("Maximum branch length: %f", *maxBrLen)
		m.SetMaxBranchLength(*maxBrLen)
		m.SetOptimizeBranchLengths()
	} else {
		log.Info("Will not optimize branch lengths")
	}

	aggMode, err := getAggModeFromString(*aggregate)
	if err != nil {
		log.Fatal(err)
	}
	if aggMode != cmodel.AggNone {
		log.Infof("Aggregation mode: %s", *aggregate)
	}
	m.SetAggregationMode(aggMode)

	if *startF != "" {
		l, err := lastLine(*startF)
		par := m.GetFloatParameters()
		if err == nil {
			err = par.ReadLine(l)
		}
		if err != nil {
			log.Debug("Reading start file as JSON")
			err2 := par.ReadFromJSON(*startF)
			// startF is neither trajectory nor correct JSON
			if err2 != nil {
				log.Error("Error reading start position from JSON:", err2)
				log.Fatal("Error reading start position from trajectory file:", err)
			}
		}
		if !par.InRange() {
			log.Fatal("Initial parameters are not in the range")
		}
	} else if *randomize {
		log.Info("Using uniform (in the boundaries) random starting point")
		par := m.GetFloatParameters()
		par.Randomize()
	}

	// iteration to skip before annealing, for adaptive mcmc
	annealingSkip := 0
	if *adaptive {
		as := optimize.NewAdaptiveSettings()
		if *skip < 0 {
			*skip = *iterations / 20
		}
		if *maxAdapt < 0 {
			*maxAdapt = *iterations / 5
		}
		annealingSkip = *maxAdapt
		log.Infof("Setting adaptive parameters, skip=%v, maxAdapt=%v", *skip, *maxAdapt)
		as.Skip = *skip
		as.MaxAdapt = *maxAdapt
		m.SetAdaptive(as)
	}

	log.Infof("Model has %d parameters.", len(m.GetFloatParameters()))

	f := os.Stdout

	if *outF != "" {
		f, err = os.Create(*outF)
		if err != nil {
			log.Fatal("Error creating trajectory file:", err)
		}
		defer f.Close()
	}

	var opt, oldOpt optimize.Optimizer

MethodLoop:
	for _, method := range strings.Split(*method, "+") {
		opt, err = getOptimizerFromString(method, *accept, annealingSkip, *seed)
		if err != nil {
			log.Error(err)
			if oldOpt == nil {
				os.Exit(1)
			} else {
				continue MethodLoop
			}
		}

		log.Infof("Using %s optimization.", method)

		opt.SetOutput(f)
		if oldOpt != nil {
			opt.LoadFromOptimizer(oldOpt)
		} else {
			opt.SetOptimizable(m)
		}
		opt.SetReportPeriod(*report)

		opt.Run(*iterations)
		summary.Optimizers = append(summary.Optimizers, opt.Summary())
		oldOpt = opt
	}
	opt.PrintFinal()
	if !*noFinal {
		m.Final()
	}
	summary.Model = m.Summary()

	if !*noOptBrLen {
		if root {
			log.Infof("unrooted_outtree=%s", t)
			err = t.Root(rootID)
			if err != nil {
				log.Error("Error rooting tree:", err)
			}
		}
		log.Infof("outtree=%s", t)
		summary.EndTree = t.String()
	}

	if *outTreeF != "" {
		f, err := os.Create(*outTreeF)
		if err != nil {
			log.Error("Error creating tree output file:", err)
		} else {
			f.WriteString(t.String() + "\n")
			f.Close()
		}
	}

	if *printFull && (aggMode != cmodel.AggNone) {
		m.SetAggregationMode(cmodel.AggNone)
		L := m.Likelihood()
		log.Notice("Full likelihood: ", L)
		summary.FullLnL = L
	}

	endTime := time.Now()

	deltaT := endTime.Sub(startTime)
	log.Noticef("Running time: %v", deltaT)
	summary.Time = deltaT.Seconds()

	// output summary in json format
	if *jsonF != "" {
		j, err := json.Marshal(summary)
		if err != nil {
			log.Error(err)
		} else {
			log.Debug(string(j))
			f, err := os.Create(*jsonF)
			if err != nil {
				log.Error("Error creating json output file:", err)
			} else {
				f.Write(j)
			}
			f.Close()
		}
	}
}

// MCMC sampler for branch-site and M0 models.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
	"syscall"
	"time"

	"github.com/op/go-logging"

	"bitbucket.org/Davydov/godon/bio"
	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/tree"
)

var githash = ""
var gitbranch = ""
var buildstamp = ""
var version = fmt.Sprintf("branch: %s, revision: %s, build time: %s", gitbranch, githash, buildstamp)

var log = logging.MustGetLogger("godon")
var formatter = logging.MustStringFormatter(`%{message}`)

func lastLine(fn string) (line string) {
	f, err := os.Open(fn)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line = scanner.Text()
	}
	if err := scanner.Err(); err != nil {
		log.Fatal("Error scanning last line:", err)
	}
	return
}

func main() {
	startTime := time.Now()

	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Godon %s\n", version)
		fmt.Fprintf(os.Stderr, "Usage:\n")
		flag.PrintDefaults()
	}

	// model
	model := flag.String("model", "M0", "todel type (M0 or BS for branch site)")
	fgBranch := flag.Int("fg", -1, "fg branch number")
	noOptBrLen := flag.Bool("nobrlen", false, "don't optimize branch lengths")
	cFreq := flag.String("cfreq", "F3X4", "codon frequecny (F0 or F3X4)")
	cFreqFileName := flag.String("cfreqfn", "", "codon frequencies file (overrides -cfreq)")
	fixw := flag.Bool("fixw", false, "fix omega=1 (for the branch-site and M8 models)")
	ncatsg := flag.Int("ncatsg", 1, "number of categories for the site gamma rate variation (no variation by default)")
	ncatcg := flag.Int("ncatcg", 1, "number of categories for the codon gamma rate variation (no variation by default)")
	ncatb := flag.Int("ncatb", 4, "number of the categories for the beta distribution (models M7&M8)")

	// optimizer parameters
	iterations := flag.Int("iter", 10000, "number of iterations")
	report := flag.Int("report", 10, "report every N iterations")
	method := flag.String("method", "simplex", "optimization method to use "+
		"(lbfgsb: limited-memory Broyden–Fletcher–Goldfarb–Shanno with bounding constraints, "+
		"simplex: downhill simplex, "+
		"annealing: simullated annealing, "+
		"mh: Metropolis-Hastings, "+
		"n_lbfgs: LBFGS from nlopt, "+
		"n_simplex: downhill simplex from nlopt, "+
		"n_cobyla: COBYLA from nlopt, "+
		"n_bobyqa: BOBYQA from nlopt, "+
		"n_sqp: SQP from nlopt"+
		"n_mlsl: MLSL from nlopt (BOBYQA local optimizer)"+
		"); you can chain optimizers with a plus sign")
	// mcmc parameters
	accept := flag.Int("accept", 200, "report acceptance rate every N iterations")

	// adaptive mcmc parameters
	adaptive := flag.Bool("adaptive", false, "use adaptive MCMC")
	skip := flag.Int("skip", -1, "number of iterations to skip for adaptive mcmc (5% by default)")
	maxAdapt := flag.Int("maxadapt", -1, "stop adapting after iteration (20% by default)")

	// optimizations
	optFixed := flag.Bool("fixed", false, "optimize fixed position likelihood by limiting markov chain to observed states")
	optAll := flag.Bool("observed", false, "optimize likelihood computation by limiting markov chain to observed states")
	printFull := flag.Bool("printfull", false, "print full likelihood in the end of optimization")

	// technical
	nCPU := flag.Int("cpu", 0, "number of cpu to use")
	seed := flag.Int64("seed", -1, "random generator seed, default time based")
	cpuProfile := flag.String("cpuprofile", "", "write cpu profile to file")

	// input/output
	outLogF := flag.String("log", "", "write log to a file")
	outF := flag.String("out", "", "write optimization trajectory to a file")
	outTreeF := flag.String("tree", "", "write tree to a file")
	startF := flag.String("start", "", "read start position from the trajectory file")
	logLevel := flag.String("loglevel", "notice", "loglevel ('critical', 'error', 'warning', 'notice', 'info', 'debug')")

	flag.Parse()

	// logging
	logging.SetFormatter(formatter)

	var backend *logging.LogBackend
	if *outLogF != "" {
		f, err := os.Create(*outLogF)
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

	// print commandline
	log.Info("Command line:", os.Args)

	if *seed == -1 {
		*seed = time.Now().UnixNano()
		log.Debug("Random seed from time")
	}
	log.Infof("Random seed=%v", *seed)

	rand.Seed(*seed)
	runtime.GOMAXPROCS(*nCPU)

	log.Infof("Using CPUs: %d.\n", runtime.GOMAXPROCS(0))

	if len(flag.Args()) < 2 {
		log.Fatal("you should specify tree and alignment")
		return
	}

	fastaFile, err := os.Open(flag.Args()[0])
	if err != nil {
		log.Fatal(err)
	}
	defer fastaFile.Close()

	ali, err := bio.ParseFasta(fastaFile)
	if err != nil {
		log.Fatal(err)
	}

	cali, err := cmodel.ToCodonSequences(ali)
	if err != nil {
		log.Fatal(err)
	}

	log.Infof("Read alignment of %d codons, %d fixed positions, %d ambiguous positions", len(cali[0].Sequence), cali.NFixed(), cali.NAmbiguous())

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
	var rootId = 0

	if t.IsRooted() {
		log.Warning("Tree is rooted. Will unroot.")
		root = true
		rootId, err = t.Unroot()
		if err != nil {
			log.Fatal("Error unrooting tree:", err)
		}
	}

	log.Infof("intree_unroot=%s", t)
	log.Debug(t.FullString())

	var cf cmodel.CodonFrequency

	if *cFreqFileName != "" {
		cFreqFile, err := os.Open(*cFreqFileName)
		if err != nil {
			log.Fatal(err)
		}
		cf, err = cmodel.ReadFrequency(cFreqFile)
		if err != nil {
			log.Fatal(err)
		}
	} else {
		switch *cFreq {
		case "F0":
			log.Info("F0 frequency")
			cf = cmodel.F0()
		case "F3X4":
			log.Info("F3X4 frequency")
			cf = cmodel.F3X4(cali)
		default:
			log.Fatal("Unknow codon freuquency specification")
		}
	}
	log.Debug(cf)

	if *fgBranch >= 0 {
		for _, node := range t.NodeIdArray() {
			if node == nil {
				continue
			}
			if node.Id == *fgBranch {
				node.Class = 1
			} else {
				node.Class = 0
			}
		}
	} else {
		class1 := 0
		for _ = range t.ClassNodes(1) {
			class1++
		}
		if class1 == 0 {
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

	var m cmodel.TreeOptimizableSiteClass

	switch *model {
	case "M0":
		log.Info("Using M0 model")
		m = cmodel.NewM0(cali, t, cf)
	case "M0vrate":
		log.Info("Using M0vrate model")
		m = cmodel.NewM0vrate(cali, t, cf)
	case "M7":
		log.Info("Using M7 model")
		log.Infof("%d beta categories, %d site gamma categories, %d codon gama categories", *ncatb, *ncatsg, *ncatcg)
		m = cmodel.NewM8(cali, t, cf, false, false, *ncatb, *ncatsg, *ncatcg)
	case "M8":
		log.Info("Using M8 model")
		log.Infof("%d beta categories, %d site gamma categories, %d codon gama categories", *ncatb, *ncatsg, *ncatcg)
		m = cmodel.NewM8(cali, t, cf, true, *fixw, *ncatb, *ncatsg, *ncatcg)
	case "BSC":
		log.Info("Using branch site C model")
		m = cmodel.NewBranchSiteC(cali, t, cf)
	case "BSG":
		log.Info("Using branch site gamma model")
		log.Infof("%d site gamma categories, %d codon gama categories", *ncatsg, *ncatcg)
		m = cmodel.NewBranchSiteGamma(cali, t, cf, *fixw, *ncatsg, *ncatcg)
	case "BS":
		log.Info("Using branch site model")
		m = cmodel.NewBranchSite(cali, t, cf, *fixw)
	default:
		log.Fatal("Unknown model specification")
	}

	log.Infof("Model has %d site class(es)", m.GetNClass())

	if !*noOptBrLen {
		log.Info("Will optimize branch lengths")
		m.SetOptimizeBranchLengths()
	} else {
		log.Info("Will not optimize branch lengths")
	}

	log.Infof("Optimize likelihood computations. Fixed positions: %t, all positions: %t.", *optFixed, *optAll)
	m.SetOptimizations(*optFixed, *optAll)

	if *startF != "" {
		l := lastLine(*startF)
		par := m.GetFloatParameters()
		err := par.ReadLine(l)
		if err != nil {
			log.Fatal("Error reading start position:", err)
		}
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
		f, err := os.Create(*outF)
		if err != nil {
			log.Fatal("Error creating trajectory file:", err)
		}
		defer f.Close()
	}

	var opt, oldOpt optimize.Optimizer

MethodLoop:
	for _, method := range strings.Split(*method, "+") {
		switch method {
		case "lbfgsb":
			lbfgsb := optimize.NewLBFGSB()
			opt = lbfgsb
		case "simplex":
			ds := optimize.NewDS()
			opt = ds
		case "mh":
			chain := optimize.NewMH(false, 0)
			chain.AccPeriod = *accept
			opt = chain
		case "annealing":
			chain := optimize.NewMH(true, annealingSkip)
			chain.AccPeriod = *accept
			opt = chain
		case "n_lbfgs":
			nlopt := optimize.NewNLOPT(optimize.NLOPT_LBFGS, *seed)
			opt = nlopt
		case "n_simplex":
			nlopt := optimize.NewNLOPT(optimize.NLOPT_SIMPLEX, *seed)
			opt = nlopt
		case "n_cobyla":
			nlopt := optimize.NewNLOPT(optimize.NLOPT_COBYLA, *seed)
			opt = nlopt
		case "n_bobyqa":
			nlopt := optimize.NewNLOPT(optimize.NLOPT_BOBYQA, *seed)
			opt = nlopt
		case "n_sqp":
			nlopt := optimize.NewNLOPT(optimize.NLOPT_SQP, *seed)
			opt = nlopt
		case "n_direct":
			nlopt := optimize.NewNLOPT(optimize.NLOPT_DIRECT, *seed)
			opt = nlopt
		case "n_crs":
			nlopt := optimize.NewNLOPT(optimize.NLOPT_CRS, *seed)
			opt = nlopt
		case "n_mlsl":
			nlopt := optimize.NewNLOPT(optimize.NLOPT_MLSL, *seed)
			opt = nlopt
		default:
			log.Errorf("Unknown optimization method: %s", method)
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
		opt.WatchSignals(os.Interrupt, syscall.SIGUSR2)

		opt.Run(*iterations)
		oldOpt = opt
	}
	opt.PrintFinal()

	if !*noOptBrLen {
		if root {
			log.Infof("unrooted_outtree=%s", t)
			err = t.Root(rootId)
			if err != nil {
				log.Error("Error rooting tree:", err)
			}
		}
		log.Infof("outtree=%s", t)
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

	if *printFull && (*optFixed || *optAll) {
		m.SetOptimizations(false, false)
		log.Notice("Full likelihood: ", m.Likelihood())
	}

	endTime := time.Now()

	log.Noticef("Running time: %v", endTime.Sub(startTime))
}

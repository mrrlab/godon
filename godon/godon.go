/*

Godon implements number of codon models (including M0 and
branch-site). It includes various likelihood optimizers as well as
Metropolis-Hastings sampler.

The basic usage of godon looks like this:

	godon alignment.fst tree.nwk

, this will run M0 model with a default optimizer (downhill simplex).

You can change a model and an optimizer:

	godon -model BS -method lbfgsb alignment.fst tree.nwk

The above will use the branch-site model and LBFGS-B optimizer.

To see all the options run:

	godon -h

*/
package main

import (
	"bufio"
	"encoding/json"
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

// Summary is storing godon summary information.
type Summary struct {
	// Version stores godon version.
	Version string `json:"version"`
	// CommandLine is an array storing binary name and all command-line parameters.
	CommandLine []string `json:"commandLine"`
	// Seed is the seed used for random number generation initialization.
	Seed int64 `json:"seed"`
	// NThreads is the number of processes used.
	NThreads int `json:"nThreads"`
	// StartTree is the starting tree after unrooting.
	StartTree string `json:"startTree"`
	// EndTree is the tree after branch length optimization (if performend).
	EndTree string `json:"endTree,omitempty"`
	// FullLnL is the full (non-aggregated) likelihood, it's only computed if specified (-printfull).
	FullLnL float64 `json:"fullLnL,omitempty"`
	// Time is the computations time in seconds.
	Time float64 `json:"time"`
	// Model is the model summary, including BEB and NEB if available.
	Model interface{} `json:"model,omitempty"`
	// Optimizers is an array of summary of all optimizers used.
	Optimizers []interface{} `json:"optimizers"`
}

// Logger settings.
var log = logging.MustGetLogger("godon")
var formatter = logging.MustStringFormatter(`%{message}`)

// lastLine returns the last line of a file content.
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

	summary := Summary{}

	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Godon %s\n", version)
		fmt.Fprintf(os.Stderr, "Usage:\n")
		flag.PrintDefaults()
	}

	// model
	model := flag.String("model", "M0", "todel type (M0 or BS for branch site)")
	gcodeId := flag.Int("gcode", 1, "NCBI genetic code id, standard by default")
	fgBranch := flag.Int("fg", -1, "fg branch number")
	maxBrLen := flag.Float64("maxbrlen", 100, "maximum branch length")
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
	startF := flag.String("start", "", "read start position from the trajectory file")
	logLevel := flag.String("loglevel", "notice", "loglevel ('critical', 'error', 'warning', 'notice', 'info', 'debug')")
	jsonF := flag.String("json", "", "write json output to a file")

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

	effectiveNCPU := runtime.GOMAXPROCS(0)
	log.Infof("Using CPUs: %d.\n", effectiveNCPU)
	summary.NThreads = effectiveNCPU

	if len(flag.Args()) < 2 {
		log.Fatal("you should specify tree and alignment")
		return
	}

	gcode, ok := bio.GeneticCodes[*gcodeId]
	if !ok {
		log.Fatalf("couldn't load genetic code with id=%d", gcodeId)
	}
	log.Infof("Genetic code: %d, \"%s\"", gcode.Id, gcode.Name)

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
	summary.StartTree = t.String()
	log.Debugf("brtree_unroot=%s", t.BrString())
	log.Debug(t.FullString())

	var cf codon.CodonFrequency

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
		if class1 == 0 &&
			(*model == "BS" || *model == "BSC" || *model == "BSG" || *model == "BSGE") {
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
	case "M1a":
		log.Info("Using M1a model")
		log.Infof("%d site gamma categories, %d codon gama categories", *ncatsg, *ncatcg)
		m = cmodel.NewM2(cali, t, cf, false, *ncatsg, *ncatcg)
	case "M2a":
		log.Info("Using M2a model")
		log.Infof("%d site gamma categories, %d codon gama categories", *ncatsg, *ncatcg)
		m = cmodel.NewM2(cali, t, cf, true, *ncatsg, *ncatcg)
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
	case "BSGE":
		log.Info("Using branch site gamma model with explicit rates")
		log.Infof("%d site gamma categories, %d codon gama categories", *ncatsg, *ncatcg)
		m = cmodel.NewBranchSiteGammaERates(cali, t, cf, *fixw, *ncatsg, *ncatcg)
	case "BS":
		log.Info("Using branch site model")
		m = cmodel.NewBranchSite(cali, t, cf, *fixw)
	default:
		log.Fatal("Unknown model specification")
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

	var aggMode cmodel.AggMode
	switch *aggregate {
	case "none":
		aggMode = cmodel.AGG_NONE
	case "observed":
		aggMode = cmodel.AGG_OBSERVED
	case "observed_new":
		aggMode = cmodel.AGG_OBSERVED_NEW
	case "fixed":
		aggMode = cmodel.AGG_FIXED
	case "random":
		aggMode = cmodel.AGG_RANDOM
	default:
		log.Fatalf("Unknown aggregation mode: %s", *aggregate)
	}
	if *aggregate != "none" {
		log.Infof("Aggregation mode: %s", *aggregate)
	}
	m.SetAggregationMode(aggMode)

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
		f, err = os.Create(*outF)
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

		opt.Run(*iterations)
		summary.Optimizers = append(summary.Optimizers, opt.Summary())
		oldOpt = opt
	}
	opt.PrintFinal()
	m.Final()
	summary.Model = m.Summary()

	if !*noOptBrLen {
		if root {
			log.Infof("unrooted_outtree=%s", t)
			err = t.Root(rootId)
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

	if *printFull && (aggMode != cmodel.AGG_NONE) {
		m.SetAggregationMode(cmodel.AGG_NONE)
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

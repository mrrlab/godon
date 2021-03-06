/*

Godon implements number of codon models (including M0 and
branch-site). It includes various likelihood optimizers as well as
Metropolis-Hastings sampler.

The basic usage of godon looks like this:

	godon M0 alignment.fst tree.nwk

, this will run M0 model with a default optimizer (LBFGS-B).

You can change a model and an optimizer:

	godon -m simplex BS alignment.fst tree.nwk

The above will use the branch-site model and downhill simplex optimizer.

To see all the options run:

	godon -h

*/
package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
	"time"

	"bitbucket.org/Davydov/godon/checkpoint"

	"gopkg.in/alecthomas/kingpin.v2"

	"github.com/op/go-logging"

	bolt "go.etcd.io/bbolt"
)

// These three variables are set during the compilation.
var githash = ""
var gitbranch = ""
var buildstamp = ""
var version = fmt.Sprintf("branch: %s, revision: %s, build time: %s", gitbranch, githash, buildstamp)

// Logger settings.
var log = logging.MustGetLogger("godon")
var formatter = logging.MustStringFormatter(`%{message}`)
var logLevels = []string{"critical", "error", "warning", "notice", "info", "debug"}

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

// warnBSGNoRate warns the user in case BSG is used with rates equal to 1.
func warnBSGNoRate() {
	if (*model == "BSG" || *hTestModel == "BSG") && *ncatsr == 1 && *ncatcr == 1 {
		log.Warning("WARNING: Attempting to use BSG model with number of categories for rate variation set to 1. This is equivalent of using BS. Perhaps you forgot to set --ncat-codon-rate.")
	}
}

// command-line options
var (
	// application
	app = kingpin.New("godon", "codon models optmizer and sampler").
		Version(version)

	opt = app.Command("optimize", "Run model optimization or sampling").Default()
	// model
	model = opt.Arg("model",
		"model type (M0 or BS for branch site)").
		Required().String()

	// optimize flags
	alignmentFileName = opt.Arg("alignment", "sequence alignment").Required().ExistingFile()
	treeFileName      = opt.Arg("tree", "starting phylogenetic tree").Required().ExistingFile()
	fixw              = opt.Flag("fix-w", "fix omega=1 (for the branch-site and M8 models)").Short('f').Bool()
	outTreeF          = opt.Flag("out-tree", "write tree to a file").String()
	printFull         = opt.Flag("full-likelihood", "print full (non-aggregated) likelihood in the end of optimization").Bool()

	// hypTest flags
	hTest      = app.Command("test", "Run test for positive selection")
	hTestModel = hTest.Arg("model",
		"model type (BS for branch site, BSG for branch-site + gamma, M2a vs M1a or M8)").
		Required().
		Enum("BS", "BSG", "M8", "M2a")
	hTestAlignmentFileName = hTest.Arg("alignment", "sequence alignment").Required().ExistingFile()
	hTestTreeFileName      = hTest.Arg("tree", "starting phylogenetic tree").Required().ExistingFile()
	sThr                   = hTest.Flag("significance-threshold",
		"LRT siginficance threshold, for H0 rerun and posterior computations").
		Default(defaultSThr).
		Float64()
	quick = hTest.Flag("quick",
		"only prevent negative LRT statistics").
		Bool()
	thorough = hTest.Flag("thorough",
		"try optimizing (e.g. reducing D-statistic) even if starting point has no improvement").
		Bool()
	m0Tree = hTest.Flag("m0-tree", "estimate branch lengths using the M0 model").
		Bool()
	m0TreeGamma = hTest.Flag("m0-tree-gamma", "use M0G model for the branch length estimation").
			Bool()
	testAllBranches = hTest.Flag("all-branches", "test all branches (for BS & BSG)").
			Bool()
	noLeavesTest = hTest.Flag("no-leaves", "don't test leaves (for BS & BSG)").
			Bool()

	//model parameters
	gcodeID       = app.Flag("gcode", "NCBI genetic code id, standard by default").Default("1").Int()
	fgBranch      = app.Flag("fg-branch", "foreground branch number").Default("-1").Int()
	maxBrLen      = app.Flag("max-branch-length", "maximum branch length").Default("100").Float64()
	noOptBrLen    = app.Flag("no-branch-length", "don't optimize branch lengths").Short('n').Bool()
	cFreq         = app.Flag("codon-frequency", "codon frequecny (F0 or F3X4)").Default("F3X4").String()
	cFreqFileName = app.Flag("codon-frequency-file", "codon frequencies file (overrides --codon-frequency)").ExistingFile()
	ncatsr        = app.Flag("ncat-site-rate", "number of categories for the site rate variation (no variation by default)").Default("1").Int()
	ncatcr        = app.Flag("ncat-codon-rate", "number of categories for the codon rate variation (no variation by default)").Default("1").Int()
	proportional  = app.Flag("proportional", "use three rates and three proportions instead of gamma distribution").Bool()
	ncatb         = app.Flag("ncat-beta", "number of the categories for the beta distribution (models M7&M8)").Default("4").Int()

	// optimizer parameters
	startF    = app.Flag("start", "read start position from the trajectory or JSON file").Short('s').ExistingFile()
	randomize = app.Flag("randomize-start", "use uniformly distributed random starting point; "+
		"by default random starting point is distributed around realistic parameter values").Bool()
	iterations = app.Flag("iter", "number of iterations").Default("10000").Int()
	report     = app.Flag("report", "report every N iterations").Default("1").Int()
	method     = app.Flag("method", "optimization method to use "+
		"(lbfgsb: limited-memory Broyden???Fletcher???Goldfarb???Shanno with bounding constraints, "+
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
		")").Short('m').Default("lbfgsb").String()

	final      = app.Flag("final", "perform final extra computations, i.e. NEB and BEB site posterior (default on, use --no-final to disable)").Default("true").Bool()
	neb        = app.Flag("neb", "perform naive empirical bayes of positive selection (default on, use --no-neb to disable)").Default("true").Bool()
	beb        = app.Flag("beb", "perform bayes empirical bayes of positive selection (default on, use --no-beb to disable)").Default("true").Bool()
	codonRates = app.Flag("codon-rates", "perform NEB analysis of codon rates").Default("false").Bool()
	siteRates  = app.Flag("site-rates", "perform NEB analysis of site rates").Default("false").Bool()
	codonOmega = app.Flag("codon-omega", "perform NEB analysis of codon omega").Default("false").Bool()

	// mcmc parameters
	accept = app.Flag("report-acceptance", "report acceptance rate every N iterations").Default("200").Int()

	// adaptive mcmc parameters
	adaptive = app.Flag("adaptive", "use adaptive MCMC or sumulated annealing").Bool()
	skip     = app.Flag("skip-adaptive", "number of iterations to skip for adaptive mcmc (5% by default)").Default("-1").Int()
	maxAdapt = app.Flag("maximum-adaptive", "stop adapting after iteration (20% by default)").Default("-1").Int()

	// optimizations
	aggregate = app.Flag("aggregate", "state aggregation mode: "+
		"observed (all positions, keep observed states), "+
		"observed_new (new implementation of observed), "+
		"fixed (absolutely conserved positions, keep observed), "+
		"random (like observed, but non-aggregated states are shuffled between the positions)").
		Default("none").Enum("none", "observed", "observed_new", "fixed", "random")
	// technical
	nThreads   = app.Flag("procs", "number of threads to use").Short('p').Int()
	seed       = app.Flag("seed", "random generator seed, default time based").Short('S').Default("-1").Int64()
	cpuProfile = app.Flag("cpu-profile", "write cpu profile to file").String()

	// input/output
	logF              = app.Flag("out", "write log to a file").Short('o').String()
	trajFn            = app.Flag("trajectory", "write optimization trajectory to a file").Short('t').String()
	checkpointFn      = app.Flag("checkpoint", "checkpoint file (possible to continue in case of termination)").Short('c').String()
	checkpointSeconds = app.Flag("checkpoint-freq", "checkpoint every N seconds (default: 30)").Default("30").Float64()
	logLevel          = app.Flag("log-level", "set loglevel "+
		"("+strings.Join(logLevels, ", ")+")").
		Short('l').Default("notice").
		Enum(logLevels...)
	jsonF = app.Flag("json", "write json output to a file").Short('j').String()

	trajF *os.File

	checkpointDB *bolt.DB
	mainBucket   = []byte("main")
)

func compareCmdLineOrSave() {
	keyCmdLine := []byte("cmdLine")
	cmdLineCp, err := checkpoint.LoadData(checkpointDB, keyCmdLine)
	if err != nil {
		log.Error("Cannot read command line from checkpoint file:", err)
	}
	cmdLineRef, err := json.Marshal(os.Args)
	if err != nil {
		log.Fatal("Cannot marshal cmdLine; this shouldn't happen:", err)
	}
	if cmdLineCp == nil {
		// command line not present in checkpoint file
		err = checkpoint.SaveData(checkpointDB, keyCmdLine, cmdLineRef)
		if err != nil {
			log.Error("Cannot save cmdLine to checkpoint file:", err)
		}
	} else {
		if !bytes.Equal(cmdLineCp, cmdLineRef) {
			log.Fatalf("Command string mismatch (checkpoint file)\n saved:   %v\n current: %v\n",
				string(cmdLineCp),
				string(cmdLineRef))
		}
	}
}

func main() {
	startTime := time.Now()
	// support -h flag
	app.HelpFlag.Short('h')
	app.VersionFlag.Short('v')
	cmd := kingpin.MustParse(app.Parse(os.Args[1:]))

	// logging
	logging.SetFormatter(formatter)

	var backend *logging.LogBackend
	if *logF != "" {
		f, err := os.OpenFile(*logF, os.O_WRONLY|os.O_CREATE|os.O_APPEND, 0666)
		if err != nil {
			log.Fatal("Error creating log file:", err)
		}
		defer f.Close()
		backend = logging.NewLogBackend(f, "", 0)
	} else {
		backend = logging.NewLogBackend(os.Stdout, "", 0)
	}
	logging.SetBackend(backend)

	level, err := logging.LogLevel(*logLevel)
	if err != nil {
		log.Fatal(err)
	}
	logging.SetLevel(level, "godon")
	logging.SetLevel(level, "optimize")
	logging.SetLevel(level, "cmodel")
	logging.SetLevel(level, "checkpoint")

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
	runtime.GOMAXPROCS(*nThreads)

	effectiveNThreads := runtime.GOMAXPROCS(0)
	log.Infof("Using threads: %d.\n", effectiveNThreads)

	if *cpuProfile != "" {
		f, err := os.Create(*cpuProfile)
		if err != nil {
			log.Fatal(err)
		}
		err = pprof.StartCPUProfile(f)
		if err != nil {
			log.Fatal(err)
		}
		defer pprof.StopCPUProfile()
	}

	callSummary := CallSummary{
		Version:     version,
		CommandLine: os.Args,
		Seed:        *seed,
		NThreads:    effectiveNThreads,
	}

	if *trajFn != "" {
		trajF, err = os.OpenFile(*trajFn, os.O_WRONLY|os.O_CREATE|os.O_APPEND, 0666)
		if err != nil {
			log.Fatalf("Error creating trajectory file: %v", err)
		}
		defer trajF.Close()
	}

	if *checkpointFn != "" {
		checkpointDB, err = bolt.Open(*checkpointFn, 0666, nil)
		defer checkpointDB.Close()
		if err != nil {
			log.Fatalf("Error opening checkpoint file: %v", err)
		}

		log.Info("Using checkpoint file")
		compareCmdLineOrSave()
	}

	var summary interface{}
	switch cmd {
	case opt.FullCommand():
		warnBSGNoRate()
		runSummary := optimization()
		summary = struct {
			*CallSummary
			OptimizationSummary
		}{&callSummary, runSummary}
	case hTest.FullCommand():
		warnBSGNoRate()
		hTestSummary, optimizations := hypTest()
		callSummary.Optimizations = optimizations
		// testAllBranches is updated by hypTest in order to
		// discriminate multiple branches vs single branch
		// tested
		if !*testAllBranches { // only a single branch was tested
			summary = struct {
				*CallSummary
				HypTestSummary
			}{&callSummary, hTestSummary[0]}
		} else { // multiple branches were tested
			callSummary.Tests = hTestSummary
			summary = &callSummary
		}
	default:
		log.Fatalf("command %v not implemented", cmd)
	}

	deltaT := time.Since(startTime)
	log.Noticef("Running time: %v", deltaT)

	callSummary.TotalTime = deltaT.Seconds()

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
				_, err := f.Write(j)
				if err != nil {
					log.Error(err)
				}
				f.Close()
			}
		}
	}
}

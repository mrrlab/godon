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
	"fmt"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
	"time"

	"gopkg.in/alecthomas/kingpin.v2"

	"github.com/op/go-logging"
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

// command-line options
var (
	// application
	app = kingpin.New("godon", "codon models optmizer and sampler").
		Version(version)

	// model
	model = app.Arg("model",
		"model type (M0 or BS for branch site)").
		Required().String()

	// input tree and alignment
	alignmentFileName = app.Arg("alignment", "sequence alignment").Required().ExistingFile()
	treeFileName      = app.Arg("tree", "starting phylogenetic tree").Required().ExistingFile()

	//model parameters
	gcodeID       = app.Flag("gcode", "NCBI genetic code id, standard by default").Default("1").Int()
	fgBranch      = app.Flag("fg-branch", "fg branch number").Default("-1").Int()
	maxBrLen      = app.Flag("max-branch-length", "maximum branch length").Default("100").Float64()
	noOptBrLen    = app.Flag("no-branch-length", "don't optimize branch lengths").Short('n').Bool()
	cFreq         = app.Flag("codon-frequency", "codon frequecny (F0 or F3X4)").Default("F3X4").String()
	cFreqFileName = app.Flag("codon-frequency-file", "codon frequencies file (overrides --codon-frequency)").ExistingFile()
	fixw          = app.Flag("fix-w", "fix omega=1 (for the branch-site and M8 models)").Short('f').Bool()
	ncatsg        = app.Flag("ncat-site-gamma", "number of categories for the site gamma rate variation (no variation by default)").Default("1").Int()
	ncatcg        = app.Flag("ncat-codon-gamma", "number of categories for the codon gamma rate variation (no variation by default)").Default("1").Int()
	ncatb         = app.Flag("ncat-beta", "number of the categories for the beta distribution (models M7&M8)").Default("4").Int()
	noFinal       = app.Flag("no-final", "don't perform final extra computations, i.e. NEB and BEB site posterior").Bool()

	// optimizer parameters
	randomize = app.Flag("randomize-start", "use uniformly distributed random starting point; "+
		"by default random starting point is distributed around realistic parameter values").Bool()
	iterations = app.Flag("iter", "number of iterations").Default("10000").Int()
	report     = app.Flag("report", "report every N iterations").Default("10").Int()
	method     = app.Flag("optimization-method", "optimization method to use "+
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
		")").Short('m').Default("lbfgsb").String()

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
	printFull = app.Flag("full-likelihood", "print full (non-aggregated) likelihood in the end of optimization").Bool()

	// technical
	nThreads   = app.Flag("procs", "number of threads to use").Short('p').Int()
	seed       = app.Flag("seed", "random generator seed, default time based").Short('S').Default("-1").Int64()
	cpuProfile = app.Flag("cpu-profile", "write cpu profile to file").String()

	// input/output
	outLogF  = app.Flag("out", "write log to a file").Short('o').String()
	outF     = app.Flag("trajectory", "write optimization trajectory to a file").Short('t').String()
	outTreeF = app.Flag("out-tree", "write tree to a file").String()
	startF   = app.Flag("start", "read start position from the trajectory or JSON file").Short('s').ExistingFile()
	logLevel = app.Flag("log-level", "set loglevel "+
		"("+strings.Join(logLevels, ", ")+")").
		Short('l').Default("notice").
		Enum(logLevels...)
	jsonF = app.Flag("json", "write json output to a file").Short('j').String()
)

func main() {
	// support -h flag
	app.HelpFlag.Short('h')
	app.VersionFlag.Short('v')
	kingpin.MustParse(app.Parse(os.Args[1:]))

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
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	summary := runOptimization(*startF, *fixw)
	summary.NThreads = effectiveNThreads
	summary.Version = version
	summary.CommandLine = os.Args
	summary.Seed = *seed

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
				f.Close()
			}
		}
	}
}

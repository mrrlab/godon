// MCMC sampler for branch-site and M0 models.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"syscall"
	"time"

	"bitbucket.org/Davydov/godon/bio"
	"bitbucket.org/Davydov/godon/cmodel"
	"bitbucket.org/Davydov/godon/optimize"
	"bitbucket.org/Davydov/godon/tree"
)

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

	// model
	model := flag.String("model", "M0", "todel type (M0 or BS for branch site)")
	fgBranch := flag.Int("fg", -1, "fg branch number")
	noOptBrLen := flag.Bool("nobrlen", false, "don't optimize branch lengths")
	cFreq := flag.String("cfreq", "F3X4", "codon frequecny (F0 or F3X4)")
	cFreqFileName := flag.String("cfreqfn", "", "codon frequencies file (overrides -cfreq)")
	fixw2 := flag.Bool("fixw2", false, "fix omega2=1 (only for branch-site model)")

	// optimizer parameters
	iterations := flag.Int("iter", 10000, "number of iterations")
	report := flag.Int("report", 10, "report every N iterations")
	method := flag.String("method", "simplex", "optimization method to use "+
		"(lbfgsb: limited-memory Broyden–Fletcher–Goldfarb–Shanno with bounding constraints, "+
		"simplex: downhill simplex, "+
		"annealing: simullated annealing, "+
		"mh: Metropolis-Hastings)")

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

	flag.Parse()

	if *outLogF != "" {
		f, err := os.Create(*outLogF)
		if err != nil {
			log.Fatal("Error creating log file:", err)
		}
		defer f.Close()
		log.SetOutput(f)
	}

	if *seed == -1 {
		*seed = time.Now().UnixNano()
		log.Println("Random seed from time")
	}
	log.Printf("Random seed=%v", *seed)

	rand.Seed(*seed)
	runtime.GOMAXPROCS(*nCPU)

	log.Printf("Using CPUs: %d.\n", runtime.GOMAXPROCS(0))

	if len(flag.Args()) < 2 {
		fmt.Println("you should specify tree and alignment")
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

	log.Printf("Read alignment of %d codons, %d fixed positions, %d ambiguous positions", len(cali[0].Sequence), cali.NFixed(), cali.NAmbiguous())

	treeFile, err := os.Open(flag.Args()[1])
	if err != nil {
		log.Fatal(err)
	}
	defer treeFile.Close()

	t, err := tree.ParseNewick(treeFile)
	if err != nil {
		log.Fatal(err)
	}

	log.Printf("intree=%s", t)
	log.Printf("brtree=%s", t.BrString())

	// Root the tree in the end
	var root = false
	var rootId = 0

	if t.IsRooted() {
		log.Print("Tree is rooted. Will unroot.")
		root = true
		rootId, err = t.Unroot()
		if err != nil {
			log.Fatal("Error unrooting tree:", err)
		}
	}

	log.Printf("intree_unroot=%s", t)
	log.Print(t.FullString())

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
			log.Print("F0 frequency")
			cf = cmodel.F0()
		case "F3X4":
			log.Print("F3X4 frequency")
			cf = cmodel.F3X4(cali)
		default:
			log.Fatal("Unknow codon freuquency specification")
		}
	}
	log.Print(cf)

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
			log.Print("Warning: no class=1 nodes")
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

	var m cmodel.TreeOptimizable

	switch *model {
	case "M0":
		log.Print("Using M0 model")
		m = cmodel.NewM0(cali, t, cf)
	case "M0vrate":
		log.Print("Using M0vrate model")
		m = cmodel.NewM0vrate(cali, t, cf)
	case "BSC":
		log.Print("Using branch site C model")
		m = cmodel.NewBranchSiteC(cali, t, cf)
	case "BSG":
		log.Print("Using branch site gamma model")
		m = cmodel.NewBranchSiteGamma(cali, t, cf, 4, *fixw2)
	case "BS":
		log.Print("Using branch site model")
		m = cmodel.NewBranchSite(cali, t, cf, *fixw2)
	default:
		log.Fatal("Unknown model specification")
	}

	if !*noOptBrLen {
		log.Print("Will optimize branch lengths")
		m.SetOptimizeBranchLengths()
	} else {
		log.Print("Will not optimize branch lengths")
	}

	log.Printf("Optimize likelihood computations. Fixed positions: %t, all positions: %t.", *optFixed, *optAll)
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
		log.Printf("Setting adaptive parameters, skip=%v, maxAdapt=%v", *skip, *maxAdapt)
		as.Skip = *skip
		as.MaxAdapt = *maxAdapt
		m.SetAdaptive(as)
	}

	log.Printf("Model has %d parameters.", len(m.GetFloatParameters()))

	var opt optimize.Optimizer
	switch *method {
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
	default:
		log.Fatal("Unknown optimization method")
	}

	log.Printf("Using %s optimization.", *method)

	opt.SetOptimizable(m)
	opt.SetReportPeriod(*report)
	opt.WatchSignals(os.Interrupt, syscall.SIGUSR2)
	if *outF != "" {
		f, err := os.Create(*outF)
		if err != nil {
			log.Fatal("Error creating trajectory file:", err)
		}
		defer f.Close()
		opt.SetOutput(f)
	}

	opt.Run(*iterations)

	if !*noOptBrLen {
		if root {
			log.Printf("unrooted_outtree=%s", t)
			err = t.Root(rootId)
			if err != nil {
				log.Print("Error rooting tree:", err)
			}
		}
		log.Printf("outtree=%s", t)
	}

	if *outTreeF != "" {
		f, err := os.Create(*outTreeF)
		if err != nil {
			log.Print("Error creating tree output file:", err)
		} else {
			f.WriteString(t.String() + "\n")
			f.Close()
		}
	}

	if *printFull && (*optFixed || *optAll) {
		m.SetOptimizations(false, false)
		log.Print("Full likelihood: ", m.Likelihood())
	}

	endTime := time.Now()

	log.Printf("Running time: %v", endTime.Sub(startTime))
}

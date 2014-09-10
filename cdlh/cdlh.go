// MCMC sampler for branch-site and M0 models.
package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"syscall"
	"time"

	"bitbucket.com/Davydov/golh/bio"
	"bitbucket.com/Davydov/golh/cmodel"
	"bitbucket.com/Davydov/golh/optimize"
	"bitbucket.com/Davydov/golh/tree"
)

func main() {
	startTime := time.Now()

	// model
	model := flag.String("model", "M0", "todel type (M0 or BS for branch site)")
	fgBranch := flag.Int("fg", -1, "fg branch number")
	noOptBrLen := flag.Bool("nobrlen", false, "don't optimize branch lengths")
	cFreq := flag.String("cfreq", "F3X4", "codon frequecny (F0 or F3X4)")
	cFreqFileName := flag.String("cfreqfn", "", "codon frequencies file (overrides -cfreq)")

	// mcmc parameters
	iterations := flag.Int("iter", 10000, "number of iterations")
	report := flag.Int("report", 10, "report every N iterations")
	accept := flag.Int("accept", 200, "report acceptance rate every N iterations")

	// adaptive mcmc parameters
	adaptive := flag.Bool("adaptive", false, "use adaptive MCMC")
	skip := flag.Int("skip", -1, "number of iterations to skip for adaptive mcmc (5% by default)")
	maxAdapt := flag.Int("maxadapt", -1, "stop adapting after iteration (20% by default)")

	// optimizations
	optFixed := flag.Bool("fixed", false, "optimize fixed position likelihood by limiting markov chain to observed states")
	optAll := flag.Bool("observed", false, "optimize likelihood computation by limiting markov chain to observed states")

	// technical
	nCPU := flag.Int("cpu", 0, "number of cpu to use")
	seed := flag.Int64("seed", -1, "random generator seed, default time based")
	cpuProfile := flag.String("cpuprofile", "", "write cpu profile to file")

	flag.Parse()

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
	treeFile, err := os.Open(flag.Args()[0])
	if err != nil {
		log.Fatal(err)
	}
	defer treeFile.Close()

	t, err := tree.ParseNewick(treeFile)
	if err != nil {
		log.Fatal(err)
	}

	log.Print(t.FullString())

	fastaFile, err := os.Open(flag.Args()[1])
	if err != nil {
		log.Fatal(err)
	}

	ali, err := bio.ParseFasta(fastaFile)
	if err != nil {
		log.Fatal(err)
	}

	cali, err := cmodel.ToCodonSequences(ali)
	if err != nil {
		log.Fatal(err)
	}

	log.Printf("Read alignment of %d codons, %d fixed positions", len(cali[0].Sequence), cali.NFixed())

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
		for _, node := range t.Nodes() {
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

	if !*noOptBrLen {
		log.Print("Will optimize branch lengths")
	} else {
		log.Print("Will not optimize branch lengths")
	}

	var m cmodel.TreeOptimizable

	switch *model {
	case "M0":
		log.Print("Using M0 model")
		m = cmodel.NewM0(cali, t, cf, !*noOptBrLen)
	case "BS":
		log.Print("Using branch site model")
		m = cmodel.NewBranchSite(cali, t, cf, !*noOptBrLen)
	default:
		log.Fatal("Unknown model specification")
	}

	log.Printf("Optimize likelihood computations. Fixed positions: %t, all positions: %t.", *optFixed, *optAll)
	m.SetOptimizations(*optFixed, *optAll)

	if *adaptive {
		as := optimize.NewAdaptiveSettings()
		if *skip < 0 {
			*skip = *iterations / 20
		}
		if *maxAdapt < 0 {
			*maxAdapt = *iterations / 5
		}
		log.Printf("Setting adaptive parameters, skip=%v, maxAdapt=%v", *skip, *maxAdapt)
		as.Skip = *skip
		as.MaxAdapt = *maxAdapt
		m.SetAdaptive(as)
	}

	log.Printf("Model has %d parameters.", len(m.GetModelParameters()))

	chain := optimize.NewMH()
	chain.AccPeriod = *accept
	var opt optimize.Optimizer = chain
	opt.SetOptimizable(m)
	opt.SetReportPeriod(*report)
	opt.WatchSignals(os.Interrupt, syscall.SIGUSR2)
	opt.Run(*iterations)

	endTime := time.Now()
	log.Printf("Runing time: %v", endTime.Sub(startTime))
}

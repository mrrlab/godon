package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	"syscall"
	"time"

	"github.com/skelterjohn/go.matrix"

	"bitbucket.com/Davydov/golh/bio"
	"bitbucket.com/Davydov/golh/mcmc"
	"bitbucket.com/Davydov/golh/tree"
)

func Sum(m *matrix.DenseMatrix) (s float64) {
	for i := 0; i < m.Rows(); i++ {
		for j := 0; j < m.Cols(); j++ {
			s += math.Abs(m.Get(i, j))
		}
	}
	return
}

func getCodons() <-chan string {
	ch := make(chan string)
	var cn func(string)
	cn = func(prefix string) {
		if len(prefix) == 3 {
			ch <- prefix
		} else {
			for _, l := range alphabet {
				cn(prefix + string(l))
			}
			if len(prefix) == 0 {
				close(ch)
			}
		}
	}
	go cn("")
	return ch
}

func codonDistance(c1, c2 string) (dist, transitions int) {
	for i := 0; i < len(c1); i++ {
		s1 := c1[i]
		s2 := c2[i]
		if s1 != s2 {
			dist++
			if ((s1 == 'A' || s1 == 'G') && (s2 == 'A' || s2 == 'G')) ||
				((s1 == 'T' || s1 == 'C') && (s2 == 'T' || s2 == 'C')) {
				transitions++
			}
		}
	}
	return
}

func createTransitionMatrix(cf CodonFrequency, kappa, omega float64, m *matrix.DenseMatrix) (*matrix.DenseMatrix, float64) {
	//fmt.Println("kappa=", kappa, ", omega=", omega)
	if m == nil {
		m = matrix.Zeros(nCodon, nCodon)
	}
	for i1 := 0; i1 < nCodon; i1++ {
		for i2 := 0; i2 < nCodon; i2++ {
			if i1 == i2 {
				m.Set(i1, i2, 0)
				continue
			}
			c1 := numCodon[byte(i1)]
			c2 := numCodon[byte(i2)]
			dist, transitions := codonDistance(c1, c2)

			if dist > 1 {
				m.Set(i1, i2, 0)
				continue
			}
			m.Set(i1, i2, cf[i2])
			if transitions == 1 {
				m.Set(i1, i2, m.Get(i1, i2)*kappa)
			}
			if bio.GeneticCode[c1] != bio.GeneticCode[c2] {
				m.Set(i1, i2, m.Get(i1, i2)*omega)
			}
		}
	}
	for i1 := 0; i1 < nCodon; i1++ {
		rowSum := 0.0
		for i2 := 0; i2 < nCodon; i2++ {
			rowSum += m.Get(i1, i2)
		}
		m.Set(i1, i1, -rowSum)
	}
	scale := 0.0
	for i := 0; i < nCodon; i++ {
		scale += -cf[i] * m.Get(i, i)
	}

	return m, scale

}

func PrintQ(Q *matrix.DenseMatrix) {
	codons := make([]string, len(codonNum))
	i := 0
	for k, _ := range codonNum {
		codons[i] = k
		i++
	}
	sort.Strings(codons)

	fmt.Print("\t")
	for _, codon := range codons {
		fmt.Print(codon, "\t")
	}
	fmt.Println()
	for _, codon1 := range codons {
		fmt.Print(codon1, "\t")
		for _, codon2 := range codons {
			fmt.Printf("%0.4f\t", Q.Get(int(codonNum[codon1]), int(codonNum[codon2])))
		}
		fmt.Println()
	}
}

func PrintUnQ(Q *matrix.DenseMatrix) {
	fmt.Print("\t")
	for i := 0; i < nCodon; i++ {
		fmt.Print(numCodon[byte(i)], "\t")
	}
	fmt.Println()
	for i1 := 0; i1 < nCodon; i1++ {
		fmt.Print(numCodon[byte(i1)], "\t")
		for i2 := 0; i2 < nCodon; i2++ {
			fmt.Printf("%0.4f\t", Q.Get(i1, i2))
		}
		fmt.Println()
	}
}

func main() {
	startTime := time.Now()
	cFreqFileName := flag.String("cfreqfn", "", "codon frequencies file (overrides -cfreq)")
	nCPU := flag.Int("cpu", 0, "number of cpu to use")
	seed := flag.Int64("seed", -1, "random generator seed, default time based")
	fgBranch := flag.Int("fg", -1, "fg branch number")
	cpuProfile := flag.String("cpuprofile", "", "write cpu profile to file")
	model := flag.String("model", "M0", "todel type (M0 or BS for branch site)")
	noOptBrLen := flag.Bool("nobrlen", false, "don't optimize branch lengths")
	cFreq := flag.String("cfreq", "F3X4", "codon frequecny (F0 or F3X4)")
	iterations := flag.Int("iter", 10000, "number of iterations")
	report := flag.Int("report", 10, "report every N iterations")
	accept := flag.Int("acceptance period", 200, "report acceptance rate every N iterations")
	adaptive := flag.Bool("adaptive", false, "use adaptive MCMC")

	flag.Parse()

	if *seed == -1 {
		*seed = time.Now().UnixNano()
		log.Println("Random seed from time")
	}
	log.Printf("Random seed=%v", *seed)

	rand.Seed(*seed)
	runtime.GOMAXPROCS(*nCPU)

	log.Printf("Using CPUs: %d.\n", runtime.GOMAXPROCS(0))
	log.Print("nCodon=", nCodon)

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

	cali, err := ToCodonSequences(ali)
	if err != nil {
		log.Fatal(err)
	}

	var cf CodonFrequency

	if *cFreqFileName != "" {
		cFreqFile, err := os.Open(*cFreqFileName)
		if err != nil {
			log.Fatal(err)
		}
		cf, err = readFrequency(cFreqFile)
		if err != nil {
			log.Fatal(err)
		}
	} else {
		switch *cFreq {
		case "F0":
			log.Print("F0 frequency")
			cf = F0()
		case "F3X4":
			log.Print("F3X4 frequency")
			cf = F3X4(cali)
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

	var m TreeOptimizable

	switch *model {
	case "M0":
		log.Print("Using M0 model")
		m = NewM0(cali, t, cf, !*noOptBrLen)
	case "BS":
		log.Print("Using branch site model")
		m = NewBranchSite(cali, t, cf, !*noOptBrLen)
	default:
		log.Fatal("Unknown model specification")
	}

	if *adaptive {
		m.SetAdaptive()
		log.Print("Setting adaptive parameters")
	}

	log.Printf("Model has %d parameters.", len(m.GetModelParameters()))

	chain := mcmc.NewMH(m)
	chain.RepPeriod = *report
	chain.AccPeriod = *accept
	chain.WatchSignals(os.Interrupt, syscall.SIGUSR2)
	chain.Run(*iterations)

	endTime := time.Now()
	log.Printf("Runing time: %v", endTime.Sub(startTime))
}

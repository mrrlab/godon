package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"

	"github.com/skelterjohn/go.matrix"

	"bitbucket.com/Davydov/golh/bio"
	"bitbucket.com/Davydov/golh/tree"
)

var (
	alphabet  = [...]byte{'T', 'C', 'A', 'G'}
	rAlphabet = map[byte]byte{'T': 0, 'C': 1, 'A': 2, 'G': 3}
	codonNum  = map[string]byte{}
	numCodon  = map[byte]string{}
	nm2id     = make(map[string]int)
	nCodon    int
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

func initCodon() int {
	i := byte(0)
	for codon := range getCodons() {
		if bio.IsStopCodon(codon) {
			continue
		}
		codonNum[codon] = i
		numCodon[i] = codon
		i++
	}
	return int(i)
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

func createTransitionMatrix(cf CodonFrequency, kappa, omega float64) (m *matrix.DenseMatrix, scale float64) {
	//fmt.Println("kappa=", kappa, ", omega=", omega)
	m = matrix.Zeros(nCodon, nCodon)
	for i1 := 0; i1 < nCodon; i1++ {
		for i2 := 0; i2 < nCodon; i2++ {
			if i1 == i2 {
				continue
			}
			c1 := numCodon[byte(i1)]
			c2 := numCodon[byte(i2)]
			dist, transitions := codonDistance(c1, c2)

			if dist > 1 {
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
	for i := 0; i < nCodon; i++ {
		scale += -cf[i] * m.Get(i, i)
	}

	return

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
	cFreqFileName := flag.String("cfreq", "", "codon frequencies file")
	nCPU := flag.Int("cpu", 0, "number of cpu to use")
	fgBranch := flag.Int("fg", -1, "fg branch number")
	cpuProfile := flag.String("cpuprofile", "", "write cpu profile to file")

	flag.Parse()

	runtime.GOMAXPROCS(*nCPU)
	mxprc := runtime.GOMAXPROCS(0)
	log.Printf("Using CPUs: %d.\n", mxprc)

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

	nCodon = initCodon()
	log.Print("nCodon=", nCodon)

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
		cf = F3X4(cali)
		//cf = F0()
	}
	log.Print(cf)

	for i, s := range cali {
		nm2id[s.Name] = i
	}

	if *fgBranch >= 0 {
		for node := range t.Nodes() {
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

	fmt.Println(M0(cali, t, cf, 2, 0.5))

	fmt.Println(H1(cali, t, cf, 2, 0.5, 0.5, 0.94702, 0.00000, 0.05298, 0.00000))
	fmt.Println(H1(cali, t, cf, 1.90991, 0.02000, 1, 0.94680, 0.00010, 0.05310, 0.00001))
	fmt.Println(H1(cali, t, cf, 1.87689, 0.01454, 1.68249, 0.89978, 0.04150, 0.05613, 0.00259))
	//fmt.Println(M3(cali, t, cf, 1.98115, 0.00424, 0.08123, 0.32679, 0.62517, 0.31211, 0.06272))
}

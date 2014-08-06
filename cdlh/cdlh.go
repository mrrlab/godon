package main

import (
	"bufio"
	"errors"
	"flag"
	"fmt"
	"math"
	"os"
	"runtime"
	"sort"
	"strconv"

	"github.com/skelterjohn/go.matrix"

	"bitbucket.com/Davydov/golh/bio"
	"bitbucket.com/Davydov/golh/tree"
)

var (
	alphabet = [...]byte{'T', 'C', 'A', 'G'}
	codonNum = map[string]byte{}
	numCodon = map[byte]string{}
	nm2id    = make(map[string]int)
	nCodon   int
)

func Sum(m *matrix.DenseMatrix) (s float64) {
	for i := 0; i < m.Rows(); i++ {
		for j := 0; j < m.Cols(); j++ {
			s += math.Abs(m.Get(i, j))
		}
	}
	return
}

func readFrequency(fileName string) (CodonFrequency, error) {
	cf := make(CodonFrequency, nCodon)
	fmt.Println("open", fileName)
	file, err := os.Open(fileName)
	if err != nil {
		return nil, err
	}

	defer file.Close()

	scanner := bufio.NewScanner(file)
	scanner.Split(bufio.ScanWords)

	codons := getCodons()
	i := 0
	for scanner.Scan() {
		codon := <-codons
		if bio.IsStopCodon(codon) {
			continue
		}
		if i >= nCodon {
			return nil, errors.New("too many frequencies in file")
		}
		f, err := strconv.ParseFloat(scanner.Text(), 64)
		if err != nil {
			return nil, err
		}
		cf[i] = f
		i++
	}
	fmt.Println(i)
	if i < nCodon {
		return nil, errors.New("not enough frequencies in file")
	}
	return cf, nil

}

func equalFrequency() CodonFrequency {
	cf := make(CodonFrequency, nCodon)
	for i := 0; i < nCodon; i++ {
		cf[i] = 1 / float64(nCodon)
	}
	return cf
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

func createTransitionMatrix(cf CodonFrequency, kappa, omega float64) (m *matrix.DenseMatrix) {
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
	scale := float64(0)
	for i := 0; i < nCodon; i++ {
		scale += -m.Get(i, i)
	}
	//fmt.Println("scale=", scale)
	m.Scale(float64(nCodon) / scale)
	scale = 0
	for i1 := 0; i1 < nCodon; i1++ {
		for i2 := 0; i2 < nCodon; i2++ {
			scale += math.Abs(m.Get(i1, i2))
		}
	}
	//fmt.Println("sum=", scale)
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

func M0(cali CodonSequences, t *tree.Tree, cf CodonFrequency, kappa, omega float64) float64 {
	Q := createTransitionMatrix(cf, kappa, omega)
	Qs := make([][]*EMatrix, t.NNodes())
	em := NewEMatrix(Q)
	err := em.Eigen()
	if err != nil {
		panic(fmt.Sprintf("error finding eigen: %v", err))
	}
	for i := 0; i < len(Qs); i++ {
		Qs[i] = make([]*EMatrix, 1)
		Qs[i][0] = em
	}
	prop := []float64{1}
	return L(cali, t, prop, Qs, cf)
}

func H1(cali CodonSequences, t *tree.Tree, cf CodonFrequency, fg int, kappa float64, omega0, omega2 float64, p0, p1, p2a, p2b float64) float64 {
	fmt.Printf("fg=%d, kappa=%f, omega0=%f, omega2=%f, p=[%f, %f, %f, %f]\n", fg, kappa, omega0, omega2, p0, p1, p2a, p2b)
	Q0 := createTransitionMatrix(cf, kappa, omega0)
	Q1 := createTransitionMatrix(cf, kappa, 1)
	Q2 := createTransitionMatrix(cf, kappa, omega2)
	em0 := NewEMatrix(Q0)
	em1 := NewEMatrix(Q1)
	em2 := NewEMatrix(Q2)
	err1 := em0.Eigen()
	err2 := em1.Eigen()
	err3 := em2.Eigen()
	if err1 != nil || err2 != nil || err3 != nil {
		panic(fmt.Sprintf("error finding eigen: %v, %v, %v", err1, err2, err3))
	}
	Qs := make([][]*EMatrix, t.NNodes())
	for i := 0; i < len(Qs); i++ {
		Qs[i] = make([]*EMatrix, 4)
		if i != fg {
			Qs[i][0] = em0
			Qs[i][1] = em1
			Qs[i][2] = em0
			Qs[i][3] = em1
		} else {
			Qs[i][0] = em0
			Qs[i][1] = em1
			Qs[i][2] = em2
			Qs[i][3] = em2
		}
	}
	return L(cali, t, []float64{p0, p1, p2a, p2b}, Qs, cf)
}

func main() {
	cFreqFile := flag.String("cfreq", "", "codon frequencies file")
	nCPU := flag.Int("cpu", 0, "number of cpu to use")

	flag.Parse()

	runtime.GOMAXPROCS(*nCPU)
	mxprc := runtime.GOMAXPROCS(0)
	fmt.Printf("Using CPUs: %d.\n", mxprc)

	if len(flag.Args()) < 2 {
		fmt.Println("you should specify tree and alignment")
		return
	}
	t, err := tree.ParseNewick(flag.Args()[0])
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Println(t.FullString())
	ali, err := bio.ParseFasta(flag.Args()[1])
	if err != nil {
		fmt.Println(err)
		return
	}

	nCodon = initCodon()
	fmt.Println("nCodon=", nCodon)

	cali, err := ToCodonSequences(ali)
	if err != nil {
		fmt.Println(err)
		return
	}

	var cf CodonFrequency

	if *cFreqFile != "" {
		cf, err = readFrequency(*cFreqFile)
		if err != nil {
			fmt.Println(err)
			return
		}
	} else {
		cf = equalFrequency()
	}
	fmt.Println(cf)

	for i, s := range cali {
		nm2id[s.Name] = i
	}

	fmt.Println(M0(cali, t, cf, 2, 0.5))

	fmt.Println(H1(cali, t, cf, 3, 2, 0.5, 0.5, 0.94702, 0.00000, 0.05298, 0.00000))
	fmt.Println(H1(cali, t, cf, 3, 1.92555, 0.02005, 1, 0.94702, 0.00000, 0.05298, 0.00000))
}

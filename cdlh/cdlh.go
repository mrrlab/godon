package main

import (
	"fmt"
	"bio"
	"os"
	"tree"
	"errors"
	"bytes"
	"bufio"
	"strconv"
	"math"
	"sort"
	"runtime"

	"github.com/skelterjohn/go.matrix"
)

var (
	alphabet = [...]byte{'T', 'C', 'A', 'G'}
	codonNum = map[string]byte{}
	numCodon = map[byte]string{}
	nm2id = make(map[string]int)	
	nCodon int
)

type CodonSequence struct {
	Name string
	Sequence []byte
}

type CodonSequences []CodonSequence

type CodonFrequency []float64

func (cf CodonFrequency) String() (s string) {
	s = "<CodonFrequency: "
	for i, f := range cf {
		s += fmt.Sprintf(" %v: %v,", numCodon[byte(i)], f)
	}
	s = s[:len(s)-1] + ">"
	return 
}

func (seq CodonSequence) String() (s string) {
	var b bytes.Buffer
	for _, c := range seq.Sequence {
		b.WriteString(numCodon[c] + " ")
	}
	s = ">" + seq.Name + "\n" + bio.Wrap(b.String(), 80)
	return
}


func (seqs CodonSequences) String() (s string) {
	for _, seq := range seqs {
		s += seq.String()
	}
	return s[:len(s)-1]
}

func Sum(m *matrix.DenseMatrix) (s float64) {
	for i := 0; i < m.Rows(); i++ {
		for j := 0; j < m.Cols(); j++ {
			s += math.Abs(m.Get(i, j))
		}
	}
	return
}

func ToCodonSequences(seqs bio.Sequences) (cs CodonSequences, err error) {
	cs = make(CodonSequences, 0, len(seqs))
	for _, seq := range seqs {
		if len(seq.Sequence) % 3 != 0 {
			return nil, errors.New("sequence length doesn't divide by 3")
		}
		var cseq CodonSequence
		cseq.Name = seq.Name
		cseq.Sequence = make([]byte, 0, len(seq.Sequence)/3)
		for i := 0; i < len(seq.Sequence); i+=3 {
			cseq.Sequence = append(cseq.Sequence, codonNum[seq.Sequence[i:i+3]])
		}
		cs = append(cs, cseq)
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
 	for ;scanner.Scan();{
		codon := <- codons
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
	for i := 0; i < nCodon; i++ {
		cf[i] = 1 / float64(nCodon)
	}
	return cf, nil

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

	/*
	go run func() {
		for _, l1 := range alphabet {
			for _, l2 := range alphabet {
				for _, l3 := range alphabet {
					ch <- string(l1) + string(l2) + string(l3)
				}
			}
		}
		close(ch)
	}()
        */
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
			dist ++
			if ((s1 == 'A' || s1 == 'G') && (s2 == 'A' || s2 == 'G')) ||
				((s1 == 'T' || s1 == 'C') && (s2 == 'T' || s2 == 'C')) {
				transitions ++
			}
		}
	}
	return
}

func createTransitionMatrix(cf CodonFrequency, kappa, omega float64) (m *matrix.DenseMatrix) {
	fmt.Println("kappa=", kappa, ", omega=", omega)
	m = matrix.Zeros(nCodon, nCodon)
	for i1 := 0; i1 < nCodon; i1++ {
		for i2 := 0; i2 < nCodon; i2++ {
			if i1 == i2 {
				continue
			}
			c1 := numCodon[byte(i1)]
			c2 := numCodon[byte(i2)]
			dist, transitions := codonDistance(c1, c2)
			
			if  dist > 1 {
				continue
			}
			m.Set(i1, i2, cf[i2])
			if transitions == 1 {
				m.Set(i1, i2, m.Get(i1, i2) * kappa)
			}
			if bio.GeneticCode[c1] != bio.GeneticCode[c2] {
				m.Set(i1, i2, m.Get(i1, i2) * omega)
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
	fmt.Println("scale=", scale)
	m.Scale(float64(nCodon) / scale)
	scale = 0
	for i1 := 0; i1 < nCodon; i1++ {
		for i2 := 0; i2 < nCodon; i2++ {
			scale += math.Abs(m.Get(i1, i2))
		}
	}
	fmt.Println("sum=", scale)
	return

}

func subL(ali CodonSequences, t *tree.Tree, eQts []*matrix.DenseMatrix, cf CodonFrequency, i int, plhch chan[][]float64) float64 {
	res := 0.0
	/*plh := make([][]float64, t.NNodes())
	for node := range t.Nodes() {
		plh[node.Id] = make([]float64, nCodon)
		plh[node.Id][0] = math.NaN()
	}*/
	plh := <- plhch
	nNodes := t.NNodes()
	if plh == nil {
		plh = make([][]float64, nNodes)
		for i := 0; i < nNodes; i ++ {
			plh[i] = make([]float64, nCodon)
		}
	}

	for i := 0; i < nNodes; i ++ {
		plh[i][0] = math.NaN()
	}

	nodes := make(chan *tree.Tree, len(ali))
	for node := range t.Terminals() {
		for l := byte(0); l < byte(nCodon); l++ {
			if l == ali[nm2id[node.Name]].Sequence[i] {
				plh[node.Id][l] = 1
			} else {
				plh[node.Id][l] = 0
			}
		}
		nodes <- node.Parent
	}

NodeLoop:
	for node := range nodes {
		for child := range node.ChildNodes() {
			if math.IsNaN(plh[child.Id][0]) {
				nodes <- node
				continue NodeLoop
			}
		}
		if ! math.IsNaN(plh[node.Id][0]) {
			continue NodeLoop
		}
		for l1 := 0; l1 < nCodon; l1 ++ {
			l := 1.0
			for child := range node.ChildNodes() {
				s := 0.0
				//fmt.Println("Q=", Q)
				//fmt.Println("Qt=", Qt)
				//fmt.Println("eQt=", eQt)
				for l2 := 0; l2 < nCodon; l2 ++ {
					s += eQts[child.Id].Get(l1, l2) * plh[child.Id][l2]
				}
				l *= s
			}
			plh[node.Id][l1] = l
		}
		nodes <- node.Parent
		if node.IsRoot() {
			for l := 0; l < nCodon; l++ { 
				res += cf[l] * plh[node.Id][l]
			}
			break NodeLoop
		}
		
	}
	close(nodes)
	plhch <- plh
	return math.Log(res)
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
	for _, codon1  := range codons{
		fmt.Print(codon1, "\t")
		for _, codon2  := range codons{
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
	for i1 := 0; i1 < nCodon; i1 ++ {
		fmt.Print(numCodon[byte(i1)], "\t")
		for i2 := 0; i2 < nCodon; i2 ++ {
			fmt.Printf("%0.4f\t", Q.Get(i1, i2))
		}
		fmt.Println()
	}
}

func L(ali CodonSequences, t *tree.Tree, Q *matrix.DenseMatrix, cf CodonFrequency) (lnL float64) {
	V, D, err := Q.Eigen()
	fmt.Println("Sum=", Sum(Q))
	if err != nil {
		panic("error finding eigen")
	}

	cD := matrix.Zeros(nCodon, nCodon)
	eQts := make([]*matrix.DenseMatrix, t.NNodes())
	iV, err := V.Inverse()
	if err != nil {
		panic("error inverting V")
	}

	for node := range t.Nodes() {
		for i := 0; i < nCodon; i ++ {
			cD.Set(i, i, math.Exp(D.Get(i, i) * node.BranchLength))
		}
		//fmt.Println(node.BranchLength, matrix.Product(V, cD, iV))
		eQts[node.Id] = matrix.Product(V, cD, iV)
	}

	// === print e^Q ===
	/*
	for i := 0; i < nCodon; i ++ {
		cD.Set(i, i, math.Exp(D.Get(i, i)))
	}
	//PrintQ(matrix.Product(V, D, iV))

	PrintUnQ(matrix.Product(V, cD, iV))

	for i1 := 0; i1 < nCodon; i1 ++ {
		for i2 := 0; i2 < nCodon; i2 ++ {
			fmt.Print(Q.Get(i1, i2), ",")
		}
	}
	fmt.Println()
	*/
	// === print e^Q ===

	mxprc := runtime.GOMAXPROCS(0)
	plhch := make(chan [][]float64, mxprc)
	fmt.Println("Using processors:", mxprc);
	for i := 0; i <  mxprc; i++ {
		plhch <- nil
	}

	ch := make(chan float64, len(ali[0].Sequence))
	fmt.Println(len(ali[0].Sequence))
	for i, _ := range ali[0].Sequence {
		go func(i int) {
			//fmt.Println(i); 
			ch <- subL(ali, t, eQts, cf, i, plhch)
		}(i);
	}


	for _, _ = range ali[0].Sequence {
		dlnL := <-ch
		//fmt.Println(dlnL)
		lnL += dlnL

		//lnL += <-ch
	}
	close(ch)
	return
}

func main() {

/*	f, err := os.Create("cpuprof")
        if err != nil {
		fmt.Println(err)
		return
        }
*/
//        pprof.StartCPUProfile(f)
//        defer pprof.StopCPUProfile()

	if (len(os.Args) < 4) {
		fmt.Println("specify files")
		return
	}
	t, err := tree.ParseNewick(os.Args[1])
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Println(t.FullString())
	ali, err := bio.ParseFasta(os.Args[2])
	if err != nil {
		fmt.Println(err)
		return
	}
	//fmt.Println(ali)

	nCodon = initCodon()
	fmt.Println("nCodon=", nCodon)

	cali, err := ToCodonSequences(ali)
	if err != nil {
		fmt.Println(err)
		return
	}
	cf, err := readFrequency(os.Args[3])
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Println(cf)

	for i, s := range cali {
		nm2id[s.Name] = i
	}

	Q := createTransitionMatrix(cf, 3, 0.5)
	fmt.Println(Q.Symmetric())
	fmt.Println(L(cali, t, Q, cf))

}


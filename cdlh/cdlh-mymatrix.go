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

	"matrix"
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
	go func() {
		for _, l1 := range alphabet {
			for _, l2 := range alphabet {
				for _, l3 := range alphabet {
					codon := string(l1) + string(l2) + string(l3)
					ch <- codon
				}
			}
		}
		close(ch)
	}()
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

func codonDistance(c1, c2 string) (dist, trans int) {
	for i := 0; i < len(c1); i++ {
		s1 := c1[i]
		s2 := c2[i]
		if s1 != s2 {
			dist ++
			if ((s1 == 'A' || s1 == 'G') && (s2 == 'A' || s2 == 'G')) ||
				((s1 == 'T' || s1 == 'C') && (s2 == 'T' || s2 == 'C')) {
				trans ++
			}
		}
	}
	return dist, trans
}

func createTransitionMatrix(cf CodonFrequency, kappa, omega float64) (m *matrix.Matrix) {
	m, _ = matrix.New(nCodon, nCodon)
	for i1 := 0; i1 < nCodon; i1++ {
		for i2 := 0; i2 < nCodon; i2++ {
			if i1 == i2 {
				continue
			}
			c1 := numCodon[byte(i1)]
			c2 := numCodon[byte(i2)]
			dist, trans := codonDistance(c1, c2)
			
			if  dist > 1 {
				continue
			}
			m.SetItem(i1, i2, cf[i2])
			if trans != 0 {
				m.ScaleItem(i1, i2, kappa)
			}
			if bio.GeneticCode[c1] != bio.GeneticCode[c2] {
				m.ScaleItem(i1, i2, omega)
			}
		}
	}
	for i1 := 0; i1 < nCodon; i1++ {
		rowSum := 0.0
		for i2 := 0; i2 < nCodon; i2++ {
			rowSum += m.GetItem(i1, i2)
		}
		m.SetItem(i1, i1, 1 - rowSum)
	}
	return

}

func isSymmetric(m *matrix.Matrix) bool {
	s1, s2 := m.GetSize()
	fmt.Println(m.GetSize)
	eps := 1E-10
	for i1 := 0; i1 < s1; i1 ++ {
		for i2 := 0; i2 < s2; i2 ++ {
			if math.Abs(m.GetItem(i1, i2) - m.GetItem(i2, i1)) > eps {
				fmt.Println(math.Abs(m.GetItem(i1, i2) - m.GetItem(i2, i1)))
				return false
			}
		}
	}
	return true
}

func subL(ali CodonSequences, t *tree.Tree, eQts map[*tree.Tree] *matrix.Matrix, cf CodonFrequency, i int) float64 {
	res := 0.0
	plh := make(map[*tree.Tree]map[byte]float64)
	for node := range t.Nodes() {
		plh[node] = make(map[byte]float64)
	}
		
	nodes := make(chan *tree.Tree, len(ali))
	for node := range t.Terminals() {
		for l := byte(0); l < byte(nCodon); l++ {
			if l == ali[nm2id[node.Name]].Sequence[i] {
				plh[node][l] = 1
			} else {
				plh[node][l] = 0
			}
		}
		nodes <- node.Parent
	}
			
			
NodeLoop:
	for node := range nodes {
		for child := range node.ChildNodes() {
			if len(plh[child]) == 0 {
				nodes <- node
				continue NodeLoop
			}
		}
		if len(plh[node]) > 0 {
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
					s += eQts[child].GetItem(l1, l2) * plh[child][byte(l2)]
				}
				l *= s
			}
			plh[node][byte(l1)] = l
		}
		nodes <- node.Parent
		if node.IsRoot() {
			for l := byte(0); l < byte(nCodon); l++ { 
				res += cf[l] * plh[node][l]
			}
			break NodeLoop
		}
		
	}
	close(nodes)
	return math.Log(res)
}


func L(ali CodonSequences, t *tree.Tree, Q *matrix.Matrix, cf CodonFrequency) (lnL float64) {
	eQts := make(map[*tree.Tree] *matrix.Matrix)
	Qt := Q.Empty()
	for node := range t.Nodes() {
		Q.Copy(Qt) // Q -> Qt
		Qt.Scale(node.BranchLength) // Qt *= t
		eQt := Q.Empty()
		Qt.Exponential(eQt)
		eQts[node] = eQt
	}

	ch := make(chan float64, len(ali[0].Sequence))
	for i, _ := range ali[0].Sequence {
		go func(i int) {
			//fmt.Println(i); 
			ch <- subL(ali, t, eQts, cf, i)
		}(i);
	}

	for _, _ = range ali[0].Sequence {
		lnL += <-ch
	}
	close(ch)
	return
}

func main() {
	if (len(os.Args) < 4) {
		fmt.Println("specify files")
		return
	}
	t, err := tree.ParseNewick(os.Args[1])
	if err != nil {
		fmt.Println(err)
		return
	}
	//fmt.Println(t.FullString())
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
	//fmt.Println(ali)
	//fmt.Println(cali)
	cf, err := readFrequency(os.Args[3])
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Println(cf)

	for i, s := range cali {
		nm2id[s.Name] = i
	}
	//fmt.Println(L(cali, t))

	Q := createTransitionMatrix(cf, 2, 0.5)
	fmt.Println(isSymmetric(Q))
	fmt.Println(L(cali, t, Q, cf))
	fmt.Println("Bye")
}


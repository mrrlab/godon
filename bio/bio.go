package bio

import (
	"fmt"
	"os"
	"strings"
	"bufio"
	"errors"
	"bytes"
)


var (
	GeneticCode = map[string]byte{
		"ATA":'I', "ATC":'I', "ATT":'I', "ATG":'M',
		"ACA":'T', "ACC":'T', "ACG":'T', "ACT":'T',
		"AAC":'N', "AAT":'N', "AAA":'K', "AAG":'K',
		"AGC":'S', "AGT":'S', "AGA":'R', "AGG":'R',
		"CTA":'L', "CTC":'L', "CTG":'L', "CTT":'L',
		"CCA":'P', "CCC":'P', "CCG":'P', "CCT":'P',
		"CAC":'H', "CAT":'H', "CAA":'Q', "CAG":'Q',
		"CGA":'R', "CGC":'R', "CGG":'R', "CGT":'R',
		"GTA":'V', "GTC":'V', "GTG":'V', "GTT":'V',
		"GCA":'A', "GCC":'A', "GCG":'A', "GCT":'A',
		"GAC":'D', "GAT":'D', "GAA":'E', "GAG":'E',
		"GGA":'G', "GGC":'G', "GGG":'G', "GGT":'G',
		"TCA":'S', "TCC":'S', "TCG":'S', "TCT":'S',
		"TTC":'F', "TTT":'F', "TTA":'L', "TTG":'L',
		"TAC":'Y', "TAT":'Y', "TAA":'_', "TAG":'_',
		"TGC":'C', "TGT":'C', "TGA":'_', "TGG":'W'}
)

func Translate(nseq string) (string, error) {
	var buffer bytes.Buffer

	if len(nseq) % 3 != 0 {
		return "", errors.New("sequence length doesn't divide by 3")
	}

	for i := 0; i < len(nseq); i+=3 {
		aa := GeneticCode[nseq[i:i+3]]
		if aa == '_' {
			if i+3 >= len(nseq) {
				// it's ok if this is the last codon
				break
			}
			return buffer.String(), errors.New("premature stop codon")
		}
		buffer.WriteByte(aa)
	}
	return buffer.String(), nil
}

func IsStopCodon(codon string) bool {
	if GeneticCode[codon] == '_' {
		return true
	}
	return false
}

type Sequence struct {
	Name string
	Sequence string
}

type Sequences []Sequence

func ParseFasta(fileName string) (seqs Sequences, err error) {
	fmt.Println("open", fileName)
	file, err := os.Open(fileName)
	if err != nil {
		return nil, err
	}

	defer file.Close()

	seqs = make(Sequences, 0, 10)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		if line[0] == '>' {
			seq := Sequence{Name: line[1:]}
			seqs = append(seqs, seq)
		} else {
			if len(seqs) == 0 {
				return nil, errors.New("sequence w/o prefix")
			}
			line = strings.ToUpper(strings.Replace(line, " ", "", -1))
			seqs[len(seqs)-1].Sequence += line
		}
	}
	return

}

func Wrap(seq string, n int) (s string) {
	for i := 0; i < len(seq); i += n {
		end := i+n
		if end > len(seq) {
			end = len(seq)
		}
		s += seq[i:end] + "\n"
	}
	return
}

func (seq Sequence) String() (s string) {
	s = ">" + seq.Name + "\n" + Wrap(seq.Sequence, 80)
	return
}


func (seqs Sequences) String() (s string) {
	for _, seq := range seqs {
		s += seq.String()
	}
	return s[:len(s)-1]
}

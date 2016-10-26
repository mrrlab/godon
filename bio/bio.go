// Package bio provides functions related to the genetic code.
package bio

import (
	"bufio"
	"bytes"
	"errors"
	"io"
	"strings"
)

var (
	// nAlphabet is a nucleotide alphabet.
	nAlphabet = [...]byte{'T', 'C', 'A', 'G'}
	// rNAlphabet is reverse nucleotide alphabet (letter to a number)
	rNAlphabet = map[byte]byte{'T': 0, 'C': 1, 'A': 2, 'G': 3}
)

// GeneticCode is a structure holding genetic code.
type GeneticCode struct {
	// Id according to NCBI.
	ID int
	// Genetic code name.
	Name string
	// Short name (if present).
	ShortName string
	// Map from codon to amino acid or '*' for stop codon.
	Map map[string]byte
	// NCodon is total number of codons (61 for the standard code).
	NCodon int
	// CodonNum translates codon into its' number.
	CodonNum map[string]byte
	// NumCodon translates codon number to its' string.
	NumCodon map[byte]string
}

// GetCodons returns a channel with every codon (64).
func GetCodons() <-chan string {
	ch := make(chan string)
	var cn func(string)
	cn = func(prefix string) {
		if len(prefix) == 3 {
			ch <- prefix
		} else {
			for _, l := range nAlphabet {
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

// newGeneticCode creates a new genetic code. id, name and shortName
// are coming from NCBI. aas is a string of amino acids, and starts is
// a string of start codons, not used currently.
func newGeneticCode(id int, name, shortName string, aas, starts string) *GeneticCode {
	gc := GeneticCode{
		ID:        id,
		Name:      name,
		ShortName: shortName,
		Map:       make(map[string]byte, 64),
		CodonNum:  make(map[string]byte, 64),
		NumCodon:  make(map[byte]string, 64),
	}
	aasBytes := []byte(aas)
	// i is the position in aas string
	i := 0
	// j is the numeric code of codons
	j := byte(0)
	for codon := range GetCodons() {
		aa := aasBytes[i]
		gc.Map[codon] = aa
		i++
		if aa != '*' {
			gc.NCodon++
			gc.CodonNum[codon] = j
			gc.NumCodon[j] = codon
			j++
		}
	}

	return &gc
}

// Translate translates nucleotide sequence string into the protein
// string. Error is returned is sequence is not divisible by three,
// non-terminal stop-codon is found or wrong codon is encountered.
func (gcode *GeneticCode) Translate(nseq string) (string, error) {
	var buffer bytes.Buffer

	if len(nseq)%3 != 0 {
		return "", errors.New("sequence length doesn't divide by 3")
	}

	// Convert all the letters to uppercase and U->T.
	nseq = strings.Replace(strings.ToUpper(nseq), "U", "T", -1)

	for i := 0; i < len(nseq); i += 3 {
		aa := gcode.Map[nseq[i:i+3]]
		if aa == 0 {
			return buffer.String(), errors.New("unknown codon")
		} else if aa == '_' {
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

// IsStopCodon tests if the string is a stop-codon (DNA alphabet,
// capital letters).
func (gcode *GeneticCode) IsStopCodon(codon string) bool {
	if gcode.Map[codon] == '*' {
		return true
	}
	return false
}

// Sequence is a type which is intended for storing nucleotide or
// protein sequence with it's name.
type Sequence struct {
	Name     string
	Sequence string
}

// Sequences stores multiple sequences. E.g. a sequence alignment.
type Sequences []Sequence

// ParseFasta parses FASTA sequences from a reader.
func ParseFasta(rd io.Reader) (seqs Sequences, err error) {
	seqs = make(Sequences, 0, 10)
	scanner := bufio.NewScanner(rd)
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

// Wrap inputs a string and wraps it so string length is n characters
// or less.
func Wrap(seq string, n int) (s string) {
	for i := 0; i < len(seq); i += n {
		end := i + n
		if end > len(seq) {
			end = len(seq)
		}
		s += seq[i:end] + "\n"
	}
	return
}

// String returns a sequence in FASTA format.
func (seq Sequence) String() (s string) {
	s = ">" + seq.Name + "\n" + Wrap(seq.Sequence, 80)
	return
}

// String returns sequences in FASTA format.
func (seqs Sequences) String() (s string) {
	for _, seq := range seqs {
		s += seq.String()
	}
	return s[:len(s)-1]
}

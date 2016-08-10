package codon

import (
	"bytes"
	"errors"
	"fmt"
	"strings"

	"bitbucket.org/Davydov/godon/bio"
)

// CodonSequence stores a sequence of codons with the sequence name.
type CodonSequence struct {
	Name     string
	Sequence []byte
}

// CodonSequences is an array (slice) of codon sequences with their
// names. E.g. codon alignment.
type CodonSequences []CodonSequence

// String returns text representation of codon frequencies.
func (cf CodonFrequency) String() (s string) {
	s = "<CodonFrequency: "
	for i, f := range cf {
		s += fmt.Sprintf(" %v: %v,", NumCodon[byte(i)], f)
	}
	s = s[:len(s)-1] + ">"
	return
}

// String returns FASTA representation of the codon sequence.
func (seq CodonSequence) String() (s string) {
	var b bytes.Buffer
	for _, c := range seq.Sequence {
		b.WriteString(NumCodon[c] + " ")
	}
	s = ">" + seq.Name + "\n" + bio.Wrap(b.String(), 80)
	return
}

// Length returns the length of codon alignment in codons.
func (seqs CodonSequences) Length() int {
	if len(seqs) == 0 {
		return 0
	}
	return len(seqs[0].Sequence)
}

// String returns all codon sequences in FASTA format.
func (seqs CodonSequences) String() (s string) {
	for _, seq := range seqs {
		s += seq.String()
	}
	return s[:len(s)-1]
}

// NAmbiguous returns number of ambiguous positions in the codon
// alignment.
func (seqs CodonSequences) NAmbiguous() (count int) {
	for i := 0; i < seqs.Length(); i++ {
		for _, seq := range seqs {
			if seq.Sequence[i] == NOCODON {
				count++
				break
			}
		}
	}
	return
}

// ToCodonSequences converts nucleotide bio.Sequences to
// CodonSequences.
func ToCodonSequences(seqs bio.Sequences) (cs CodonSequences, err error) {
	cs = make(CodonSequences, 0, len(seqs))
	for _, seq := range seqs {
		if len(seq.Sequence)%3 != 0 {
			return nil, errors.New("sequence length doesn't divide by 3")
		}
		var cseq CodonSequence
		cseq.Name = seq.Name
		cseq.Sequence = make([]byte, 0, len(seq.Sequence)/3)
		for i := 0; i < len(seq.Sequence); i += 3 {
			cnum, ok := CodonNum[strings.ToUpper(seq.Sequence[i:i+3])]
			if !ok {
				cnum = NOCODON
			}
			cseq.Sequence = append(cseq.Sequence, cnum)
		}
		cs = append(cs, cseq)
	}
	return
}

// NFixed calculates number of constant positions in the alignment.
func (seqs CodonSequences) NFixed() (f int) {
	f = seqs.Length()
	for pos := 0; pos < seqs.Length(); pos++ {
		for i := 1; i < len(seqs); i++ {
			if seqs[i].Sequence[pos] != seqs[0].Sequence[pos] {
				f--
				break
			}
		}
	}
	return
}

// Fixed returns a bool vector, all absolutely conserved (fixed)
// positions have true value.
func (seqs CodonSequences) Fixed() (fixed []bool) {
	fixed = make([]bool, seqs.Length())
	for pos := 0; pos < seqs.Length(); pos++ {
		isFixed := true
		for i := 1; i < len(seqs); i++ {
			if seqs[i].Sequence[pos] != seqs[0].Sequence[pos] {
				isFixed = false
				break
			}
		}
		if isFixed {
			fixed[pos] = true
		}
	}
	return
}

// Letters returns a set of present and absent codons at each position
// of the alignment.
func (seqs CodonSequences) Letters() (found [][]int, absent [][]int) {
	found = make([][]int, seqs.Length())
	absent = make([][]int, seqs.Length())
	for pos := 0; pos < seqs.Length(); pos++ {
		found[pos] = make([]int, 0, NCodon)
		absent[pos] = make([]int, 0, NCodon)
		pf := make(map[int]bool, NCodon)
		for i := 0; i < len(seqs); i++ {
			if seqs[i].Sequence[pos] != NOCODON {
				pf[int(seqs[i].Sequence[pos])] = true
			}
		}
		for l := 0; l < NCodon; l++ {
			if pf[l] {
				found[pos] = append(found[pos], l)
			} else {
				absent[pos] = append(absent[pos], l)
			}
		}

		if len(found[pos]) < NCodon {
			found[pos] = append(found[pos], NCodon)
		}
	}
	return
}

// ExtendedLetters returns a set of codons for all present aminoacids
// and absent ones at each position of the alignment.
func (seqs CodonSequences) ExtendedLetters() (found [][]int, absent [][]int) {
	found = make([][]int, seqs.Length())
	absent = make([][]int, seqs.Length())
	for pos := 0; pos < seqs.Length(); pos++ {
		found[pos] = make([]int, 0, NCodon)
		absent[pos] = make([]int, 0, NCodon)
		pf := make(map[int]bool, NCodon)
		for i := 0; i < len(seqs); i++ {
			if seqs[i].Sequence[pos] != NOCODON {
				codon := seqs[i].Sequence[pos]
				scodon := NumCodon[codon]
				aa := bio.GeneticCode[scodon]
				for _, c := range bio.RGeneticCode[aa] {
					pf[int(CodonNum[c])] = true
				}
			}
		}
		for l := 0; l < NCodon; l++ {
			if pf[l] {
				found[pos] = append(found[pos], l)
			} else {
				absent[pos] = append(absent[pos], l)
			}
		}

		if len(found[pos]) < NCodon {
			found[pos] = append(found[pos], NCodon)
		}
	}
	return
}

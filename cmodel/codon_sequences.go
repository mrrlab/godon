package cmodel

import (
	"bytes"
	"errors"
	"fmt"

	"bitbucket.com/Davydov/godon/bio"
)

var (
	alphabet  = [...]byte{'T', 'C', 'A', 'G'}
	rAlphabet = map[byte]byte{'T': 0, 'C': 1, 'A': 2, 'G': 3}
	codonNum  = map[string]byte{}
	numCodon  = map[byte]string{}
	nCodon    int
)

const (
	NOCODON = 255
)

type CodonSequence struct {
	Name     string
	Sequence []byte
}

type CodonSequences []CodonSequence

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

func (seqs CodonSequences) NAmbiguous() (count int) {
	for i := range seqs[0].Sequence {
		for _, seq := range seqs {
			if seq.Sequence[i] == NOCODON {
				count++
				break
			}
		}
	}
	return
}

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
			cnum, ok := codonNum[seq.Sequence[i:i+3]]
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
	f = len(seqs[0].Sequence)
	for pos, _ := range seqs[0].Sequence {
		for i := 1; i < len(seqs); i++ {
			if seqs[i].Sequence[pos] != seqs[0].Sequence[pos] {
				f--
				break
			}
		}
	}
	return
}

func (seqs CodonSequences) Fixed() (fixed []bool) {
	fixed = make([]bool, len(seqs[0].Sequence))
	for pos, _ := range seqs[0].Sequence {
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

func (seqs CodonSequences) Letters() (found [][]int, absent [][]int) {
	found = make([][]int, len(seqs[0].Sequence))
	absent = make([][]int, len(seqs[0].Sequence))
	for pos, _ := range seqs[0].Sequence {
		found[pos] = make([]int, 0, nCodon)
		absent[pos] = make([]int, 0, nCodon)
		pf := make(map[int]bool, nCodon)
		for i := 0; i < len(seqs); i++ {
			pf[int(seqs[i].Sequence[pos])] = true
		}
		for l := 0; l < nCodon; l++ {
			if pf[l] {
				found[pos] = append(found[pos], l)
			} else {
				absent[pos] = append(absent[pos], l)
				//found[pos] = append(found[pos], l)
			}
		}

		if len(found[pos]) < nCodon {
			found[pos] = append(found[pos], nCodon)
		}
	}
	return
}

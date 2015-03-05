package cmodel

import (
	"bufio"
	"errors"
	"io"
	"strconv"

	"bitbucket.com/Davydov/godon/bio"
)

func init() {
	i := byte(0)
	for codon := range getCodons() {
		if bio.IsStopCodon(codon) {
			continue
		}
		codonNum[codon] = i
		numCodon[i] = codon
		i++
	}
	nCodon = int(i)
}

type CodonFrequency []float64

// getCodons returns a channel with every codon (64).
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

func ReadFrequency(rd io.Reader) (CodonFrequency, error) {
	cf := make(CodonFrequency, nCodon)

	scanner := bufio.NewScanner(rd)
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
	if i < nCodon {
		return nil, errors.New("not enough frequencies in file")
	}
	return cf, nil

}

func F0() CodonFrequency {
	cf := make(CodonFrequency, nCodon)
	for i := 0; i < nCodon; i++ {
		cf[i] = 1 / float64(nCodon)
	}
	return cf
}

func F3X4(cali CodonSequences) (cf CodonFrequency) {
	poscf := make([][]float64, 3)
	for i := 0; i < 3; i++ {
		poscf[i] = make([]float64, 4)
	}

	for _, cs := range cali {
		for _, codon := range cs.Sequence {
			if codon == NOCODON {
				continue
			}
			cs := numCodon[byte(codon)]
			poscf[0][rAlphabet[cs[0]]]++
			poscf[1][rAlphabet[cs[1]]]++
			poscf[2][rAlphabet[cs[2]]]++
		}
	}

	cf = make(CodonFrequency, nCodon)

	sum := float64(0)
	for ci, cs := range numCodon {
		cf[ci] = poscf[0][rAlphabet[cs[0]]] * poscf[1][rAlphabet[cs[1]]] * poscf[2][rAlphabet[cs[2]]]
		sum += cf[ci]
	}

	for ci, _ := range cf {
		cf[ci] /= sum
	}

	return
}

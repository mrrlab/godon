package main

import (
	"bufio"
	"errors"
	"io"
	"strconv"

	"bitbucket.com/Davydov/golh/bio"
)

type CodonFrequency []float64

func readFrequency(rd io.Reader) (CodonFrequency, error) {
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

func equalFrequency() CodonFrequency {
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

package main

import (
	"io"
	"bufio"
	"strconv"
	"errors"

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

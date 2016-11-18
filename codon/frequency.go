package codon

import (
	"bufio"
	"errors"
	"io"
	"strconv"

	"bitbucket.org/Davydov/godon/bio"
)

var (
	// rAlphabet is reverse nucleotide alphabet (letter to a number)
	rAlphabet = map[byte]byte{'T': 0, 'C': 1, 'A': 2, 'G': 3}
)

// Frequency is array (slice) of codon frequencies.
type Frequency struct {
	Freq  []float64
	GCode *bio.GeneticCode
}

// ReadFrequency reads codon frequencies from a reader. It should be
// just a list of numbers in a text format.
func ReadFrequency(rd io.Reader, gcode *bio.GeneticCode) (Frequency, error) {
	cf := Frequency{
		Freq:  make([]float64, gcode.NCodon),
		GCode: gcode,
	}

	scanner := bufio.NewScanner(rd)
	scanner.Split(bufio.ScanWords)

	codons := bio.GetCodons()
	i := 0
	for scanner.Scan() {
		codon := <-codons
		if gcode.IsStopCodon(codon) {
			continue
		}
		if i >= gcode.NCodon {
			return cf, errors.New("too many frequencies in file")
		}
		f, err := strconv.ParseFloat(scanner.Text(), 64)
		if err != nil {
			return cf, err
		}
		cf.Freq[i] = f
		i++
	}
	if i < gcode.NCodon {
		return cf, errors.New("not enough frequencies in file")
	}
	return cf, nil

}

// F0 returns array (slice) of equal codon frequencies.
func F0(cali Sequences) Frequency {
	gcode := cali[0].GCode
	cf := Frequency{
		Freq:  make([]float64, gcode.NCodon),
		GCode: gcode,
	}
	for i := 0; i < gcode.NCodon; i++ {
		cf.Freq[i] = 1 / float64(gcode.NCodon)
	}
	return cf
}

// F3X4 computes F3X4-style frequencies based on the alignment.
func F3X4(cali Sequences) (cf Frequency) {
	gcode := cali[0].GCode
	poscf := make([][]float64, 3)
	for i := 0; i < 3; i++ {
		poscf[i] = make([]float64, 4)
	}

	for _, cs := range cali {
		for _, codon := range cs.Sequence {
			if codon == NOCODON {
				continue
			}
			cs := gcode.NumCodon[codon]
			poscf[0][rAlphabet[cs[0]]]++
			poscf[1][rAlphabet[cs[1]]]++
			poscf[2][rAlphabet[cs[2]]]++
		}
	}

	cf = Frequency{
		Freq:  make([]float64, gcode.NCodon),
		GCode: gcode,
	}

	sum := float64(0)
	for ci, cs := range gcode.NumCodon {
		cf.Freq[ci] = poscf[0][rAlphabet[cs[0]]] * poscf[1][rAlphabet[cs[1]]] * poscf[2][rAlphabet[cs[2]]]
		sum += cf.Freq[ci]
	}

	for ci := range cf.Freq {
		cf.Freq[ci] /= sum
	}

	return
}

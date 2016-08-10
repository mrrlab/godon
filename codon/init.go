// Codon is the package for working with codons, codon transition
// matrices and codon frequencies.
package codon

import (
	"github.com/skelterjohn/go.matrix"

	"bitbucket.org/Davydov/godon/bio"
)

var (
	// alphabet is a nucleotide alphabet.
	alphabet = [...]byte{'T', 'C', 'A', 'G'}
	// rAlphabet is reverse nucleotide alphabet (letter to a number)
	rAlphabet = map[byte]byte{'T': 0, 'C': 1, 'A': 2, 'G': 3}
	// CodonNum translates codon into its' number.
	CodonNum = map[string]byte{}
	// NumCodon translates codon number to its' string.
	NumCodon = map[byte]string{}
	// NCodon is total number of codons (61).
	NCodon int
	// NOCODON is a special number which represents an unknown
	// codon.
	NOCODON = byte(255)
)

// init initializes codon to a number arrays, number of codons and
// zero Q and identity P matrices.
func init() {
	i := byte(0)
	for codon := range getCodons() {
		if bio.IsStopCodon(codon) {
			continue
		}
		CodonNum[codon] = i
		NumCodon[i] = codon
		i++
	}
	NCodon = int(i)
	ZeroQ = matrix.Zeros(NCodon, NCodon)
	IdentityP = createIdentityMatrix(NCodon)
}

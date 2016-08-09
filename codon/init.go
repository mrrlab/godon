// Codon is the package for working with codons, codon transition
// matrices and codon frequencies.
package codon

import (
	"github.com/skelterjohn/go.matrix"

	"bitbucket.org/Davydov/godon/bio"
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

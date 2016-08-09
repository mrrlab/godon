// Codon is the package for working with codons, codon transition
// matrices and codon frequencies.
package codon

import (
	"github.com/gonum/matrix/mat64"

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
	ZeroQ = mat64.NewDense(NCodon, NCodon, nil)
	IdentityP = createIdentityMatrix(NCodon)
}

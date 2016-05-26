package codon

import (
	"github.com/skelterjohn/go.matrix"

	"bitbucket.org/Davydov/godon/bio"
)

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

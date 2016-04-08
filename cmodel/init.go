package cmodel

import (
	"github.com/gonum/matrix/mat64"

	"bitbucket.org/Davydov/godon/bio"
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
	zeroQ = mat64.NewDense(nCodon, nCodon, nil)
	identityP = createIdentityMatrix(nCodon)

}

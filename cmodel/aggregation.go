package cmodel

import (
	"math"

	"bitbucket.org/Davydov/godon/codon"
)

// AggMode is a type specifying aggregation mode.
type AggMode int

// Aggregation modes.
const (
	// No aggregation.
	AGG_NONE AggMode = iota
	// Aggregation on all the positions. All non-observed states
	// are aggregated.
	AGG_OBSERVED
	// Aggregation on absolutely conserved positions. All
	// non-observed states are aggregated.
	AGG_FIXED
	// Aggregation on all the positions. All non-observed states
	// are aggregated. More general implementation.
	AGG_OBSERVED_NEW
	// Aggregation on all the positions. Similar to observed, but
	// a set of non-aggregated states is shuffled between the
	// alignment positions.
	AGG_RANDOM
)

// observedSubL calculates likelihood for given site class and position
// taking into account only visible states.
func (m *BaseModel) observedSubL(class, pos int, plh [][]float64, lettersF, lettersA []int) (res float64) {
	if len(lettersA) <= 1 {
		// aggregation makes sense only for two absent
		// letters or more
		return m.fullSubL(class, pos, plh)
	}
	fabs := 0.0
	for _, l := range lettersA {
		fabs += m.cf[l]
	}

	for node := range m.tree.Terminals() {
		cod := m.cali[node.LeafId].Sequence[pos]
		nfound := 0
		for _, l := range lettersF {
			if l == int(m.cali[node.LeafId].Sequence[pos]) || cod == codon.NOCODON {
				plh[node.Id][l] = 1
				nfound++
			} else {
				plh[node.Id][l] = 0
			}
		}
		if cod == codon.NOCODON || nfound == 0 {
			plh[node.Id][codon.NCodon] = 1
		}
	}

	for _, node := range m.tree.NodeOrder() {
		for _, l1 := range lettersF {
			l := 1.0
			for _, child := range node.ChildNodes() {
				// get child partial likelhiood
				cplh := plh[child.Id]
				s := 0.0
				if l1 != codon.NCodon {
					// get the row
					q := m.eQts[class][child.Id][l1*codon.NCodon:]

					for _, l2 := range lettersF {
						//s += q.Get(l1, l2) * plh[child.Id][l2]
						if l2 != codon.NCodon {
							s += q[l2] * cplh[l2]
						} else {
							pia := 0.0
							for _, l2 := range lettersA {
								pia += q[l2]
							}
							s += pia * cplh[l2]
						}
					}
				} else {

					paa := 1.0
					for _, l2 := range lettersF {
						pai := 0.0
						if l2 != codon.NCodon {
							for _, l1 := range lettersA {
								pai += m.cf[l1] * m.eQts[class][child.Id][l1*codon.NCodon+l2]
							}
							pai /= fabs

							paa -= pai
							s += pai * cplh[l2]
						} else {
							s += paa * cplh[l2]
						}
					}
				}
				l *= s
			}
			plh[node.Id][l1] = l
		}

		if node.IsRoot() {
			for _, l := range lettersF {
				if l != codon.NCodon {

					res += m.cf[l] * plh[node.Id][l]
				} else {
					res += fabs * plh[node.Id][l]
				}
			}
			break
		}

	}
	return
}

// observedStates creates codon2state, state2codons translation slices
// & state frequencies for observed codon-based aggregation.
func (m *BaseModel) observedStates(lettersF, lettersA []int) (codon2state []int, state2codons [][]int, stateFreq []float64) {
	NStates := len(lettersF)

	codon2state = make([]int, codon.NCodon)
	state2codons = make([][]int, NStates)
	stateFreq = make([]float64, NStates)
	for i, l := range lettersF {
		if l != codon.NCodon {
			state2codons[i] = append(state2codons[i], l)
			codon2state[l] = i
			stateFreq[i] += m.cf[l]
		}
	}
	aState := NStates - 1
	for _, l := range lettersA {
		state2codons[aState] = append(state2codons[aState], l)
		codon2state[l] = aState
		stateFreq[aState] += m.cf[l]
	}
	return
}

// aggSubL calculates likelihood for given site class and position
// using provided aggregation scheme.
func (m *BaseModel) aggSubL(class, pos int, plh [][]float64, codon2state []int, state2codons [][]int, stateFreq []float64) (res float64) {
	NStates := len(state2codons)
	for node := range m.tree.Terminals() {
		cod := m.cali[node.LeafId].Sequence[pos]
		st := codon2state[cod]
		for i := range state2codons {
			if cod == codon.NOCODON || i == st {
				plh[node.Id][i] = 1
			} else {
				plh[node.Id][i] = 0
			}
		}
	}

	for _, node := range m.tree.NodeOrder() {
		for s1 := 0; s1 < NStates; s1++ {
			l := 1.0
			for _, child := range node.ChildNodes() {
				// get child partial likelhiood
				cplh := plh[child.Id]
				s := 0.0
				for s2 := 0; s2 < NStates; s2++ {
					ps12 := 0.0
					for _, l1 := range state2codons[s1] {
						// get the row
						q := m.eQts[class][child.Id][l1*codon.NCodon:]
						pl12 := 0.0
						for _, l2 := range state2codons[s2] {
							pl12 += q[l2]

						}
						ps12 += m.cf[l1] * pl12
					}

					//s += q.Get(l1, l2) * plh[child.Id][l2]
					s += ps12 / stateFreq[s1] * cplh[s2]
				}
				l *= s
			}
			plh[node.Id][s1] = l
		}

		if node.IsRoot() {
			for s := 0; s < NStates; s++ {
				res += stateFreq[s] * plh[node.Id][s]
			}
			break
		}

	}

	return
}

// fixedSubL calculates likelihood for given site class and position
// if the site is fixed.
func (m *BaseModel) fixedSubL(class, pos int, plh [][]float64) (res float64) {
	for i := 0; i < m.tree.MaxNodeId()+1; i++ {
		plh[i][0] = math.NaN()
	}

	l := m.lettersF[pos][0]

	for node := range m.tree.Terminals() {
		plh[node.Id][0] = 1
		plh[node.Id][1] = 0
	}

	for _, node := range m.tree.NodeOrder() {
		plh[node.Id][0] = 1
		plh[node.Id][1] = 1
		for _, child := range node.ChildNodes() {
			p00 := m.eQts[class][child.Id][l*codon.NCodon+l]
			p01 := 1 - p00
			p10 := 0.0
			for l1 := 0; l1 < codon.NCodon; l1++ {
				if l != l1 {
					p10 += m.cf[l1] * m.eQts[class][child.Id][l1*codon.NCodon+l]
				}
			}
			p10 /= (1 - m.cf[l])
			p11 := 1 - p10

			cplh := plh[child.Id]
			plh[node.Id][0] *= p00*cplh[0] + p01*cplh[1]
			plh[node.Id][1] *= p10*cplh[0] + p11*cplh[1]
		}

		if node.IsRoot() {
			res = m.cf[l]*plh[node.Id][0] + (1-m.cf[l])*plh[node.Id][1]
			break
		}

	}
	return
}

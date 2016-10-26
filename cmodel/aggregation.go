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
	AggNone AggMode = iota
	// Aggregation on all the positions. All non-observed states
	// are aggregated.
	AggObserved
	// Aggregation on absolutely conserved positions. All
	// non-observed states are aggregated.
	AggFixed
	// Aggregation on all the positions. All non-observed states
	// are aggregated. More general implementation.
	AggObservedNew
	// Aggregation on all the positions. Similar to observed, but
	// a set of non-aggregated states is shuffled between the
	// alignment positions.
	AggRandom
)

// aggSchema defines an aggregation schema.
type aggSchema struct {
	codon2state  []int
	state2codons [][]int
	stateFreq    []float64
}

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
		fabs += m.cf.Freq[l]
	}

	NCodon := m.cf.GCode.NCodon

	for node := range m.tree.Terminals() {
		cod := m.cali[node.LeafID].Sequence[pos]
		nfound := 0
		for _, l := range lettersF {
			if l == int(m.cali[node.LeafID].Sequence[pos]) || cod == codon.NOCODON {
				plh[node.ID][l] = 1
				nfound++
			} else {
				plh[node.ID][l] = 0
			}
		}
		if cod == codon.NOCODON || nfound == 0 {
			plh[node.ID][NCodon] = 1
		}
	}

	for _, node := range m.tree.NodeOrder() {
		for _, l1 := range lettersF {
			l := 1.0
			for _, child := range node.ChildNodes() {
				// get child partial likelhiood
				cplh := plh[child.ID]
				s := 0.0
				if l1 != NCodon {
					// get the row
					q := m.eQts[class][child.ID][l1*NCodon:]

					for _, l2 := range lettersF {
						//s += q.Get(l1, l2) * plh[child.Id][l2]
						if l2 != NCodon {
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
						if l2 != NCodon {
							for _, l1 := range lettersA {
								pai += m.cf.Freq[l1] * m.eQts[class][child.ID][l1*NCodon+l2]
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
			plh[node.ID][l1] = l
		}

		if node.IsRoot() {
			for _, l := range lettersF {
				if l != NCodon {

					res += m.cf.Freq[l] * plh[node.ID][l]
				} else {
					res += fabs * plh[node.ID][l]
				}
			}
			break
		}

	}
	return
}

// observedStates creates codon2state, state2codons translation slices
// & state frequencies for observed codon-based aggregation.
func (m *BaseModel) observedStates(lettersF, lettersA []int) (schema *aggSchema) {
	NStates := len(lettersF)

	schema = &aggSchema{
		codon2state:  make([]int, m.cf.GCode.NCodon),
		state2codons: make([][]int, NStates),
		stateFreq:    make([]float64, NStates),
	}

	for i, l := range lettersF {
		if l != m.cf.GCode.NCodon {
			schema.state2codons[i] = append(schema.state2codons[i], l)
			schema.codon2state[l] = i
			schema.stateFreq[i] += m.cf.Freq[l]
		}
	}
	aState := NStates - 1
	for _, l := range lettersA {
		schema.state2codons[aState] = append(schema.state2codons[aState], l)
		schema.codon2state[l] = aState
		schema.stateFreq[aState] += m.cf.Freq[l]
	}
	return
}

// aggSubL calculates likelihood for given site class and position
// using provided aggregation schema.
func (m *BaseModel) aggSubL(class, pos int, plh [][]float64, schema *aggSchema) (res float64) {
	NStates := len(schema.state2codons)
	NCodon := m.cf.GCode.NCodon

	for node := range m.tree.Terminals() {
		cod := m.cali[node.LeafID].Sequence[pos]
		st := schema.codon2state[cod]
		for i := range schema.state2codons {
			if cod == codon.NOCODON || i == st {
				plh[node.ID][i] = 1
			} else {
				plh[node.ID][i] = 0
			}
		}
	}

	for _, node := range m.tree.NodeOrder() {
		for s1 := 0; s1 < NStates; s1++ {
			l := 1.0
			for _, child := range node.ChildNodes() {
				// get child partial likelhiood
				cplh := plh[child.ID]
				s := 0.0
				rowSum := 0.0
				for s2 := 0; s2 < NStates; s2++ {
					// probability of substitution s1 -> s2
					ps12 := 0.0
					if s2 != NStates-1 {
						for _, l1 := range schema.state2codons[s1] {
							// get the row
							q := m.eQts[class][child.ID][l1*NCodon:]
							pl12 := 0.0
							for _, l2 := range schema.state2codons[s2] {
								pl12 += q[l2]

							}
							ps12 += m.cf.Freq[l1] * pl12
						}
						//s += q.Get(l1, l2) * plh[child.Id][l2]
						ps12 /= schema.stateFreq[s1]
						rowSum += ps12
					} else {
						// sum of probabilities should be equal to 1.
						// we do not have to compute the full probability.
						// this saves a lot of time, especially if the
						// larst states includes many codons.
						ps12 = 1 - rowSum
					}
					s += ps12 * cplh[s2]
				}
				l *= s
			}
			plh[node.ID][s1] = l
		}

		if node.IsRoot() {
			for s := 0; s < NStates; s++ {
				res += schema.stateFreq[s] * plh[node.ID][s]
			}
			break
		}

	}

	return
}

// fixedSubL calculates likelihood for given site class and position
// if the site is fixed.
func (m *BaseModel) fixedSubL(class, pos int, plh [][]float64) (res float64) {
	NCodon := m.cf.GCode.NCodon

	for i := 0; i < m.tree.MaxNodeID()+1; i++ {
		plh[i][0] = math.NaN()
	}

	l := m.lettersF[pos][0]

	for node := range m.tree.Terminals() {
		plh[node.ID][0] = 1
		plh[node.ID][1] = 0
	}

	for _, node := range m.tree.NodeOrder() {
		plh[node.ID][0] = 1
		plh[node.ID][1] = 1
		for _, child := range node.ChildNodes() {
			p00 := m.eQts[class][child.ID][l*NCodon+l]
			p01 := 1 - p00
			p10 := 0.0
			for l1 := 0; l1 < NCodon; l1++ {
				if l != l1 {
					p10 += m.cf.Freq[l1] * m.eQts[class][child.ID][l1*NCodon+l]
				}
			}
			p10 /= (1 - m.cf.Freq[l])
			p11 := 1 - p10

			cplh := plh[child.ID]
			plh[node.ID][0] *= p00*cplh[0] + p01*cplh[1]
			plh[node.ID][1] *= p10*cplh[0] + p11*cplh[1]
		}

		if node.IsRoot() {
			res = m.cf.Freq[l]*plh[node.ID][0] + (1-m.cf.Freq[l])*plh[node.ID][1]
			break
		}

	}
	return
}

package main

// copyMap creates a copy of the map
func copyMap(m map[string]float64) (r map[string]float64) {
	r = make(map[string]float64, len(m))
	for k, v := range m {
		r[k] = v
	}
	return
}

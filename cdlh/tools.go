package main

func Reflect(val, min, max float64) float64 {
	if min >= max {
		return val
	}
	for val < min || val > max {
		if val < min {
			val = min + (min - val)
		}
		if val > max {
			val = max - (val - max)
		}
	}
	return val
}

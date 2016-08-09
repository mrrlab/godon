package optimize

import (
	"bufio"
	"strconv"
	"strings"
)

// ReadFloats converts string of floats into slice of float64.
func ReadFloats(s string) ([]float64, error) {
	r := strings.NewReader(s)
	scanner := bufio.NewScanner(r)
	scanner.Split(bufio.ScanWords)
	var result []float64
	for scanner.Scan() {
		x, err := strconv.ParseFloat(scanner.Text(), 64)
		if err != nil {
			return result, err
		}
		result = append(result, x)
	}
	return result, scanner.Err()
}

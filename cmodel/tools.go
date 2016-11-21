package cmodel

import (
	"bytes"
	"fmt"
	"path"

	"github.com/op/go-logging"
)

const (
	data1 = "EMGT00050000008747.Drosophila.002"
	data2 = "ENSGT00550000073950.Euteleostomi.07.001"
	data3 = "EMGT00050000000025.Drosophila.001"

	// smallDiff is a threshold for testing
	// if likelihood ratio is larger error is emmited
	smallDiff = 1e-3
)

// log is a global logging variable.
var log = logging.MustGetLogger("cmodel")

// GetTreeAlignment returns a tree and alignment for testing purposes.
func GetTreeAlignment(dataset string, cfreq string) (*Data, error) {
	alifn := path.Join("testdata", dataset+".fst")
	treefn := path.Join("testdata", dataset+".nwk")

	// using standard genetic code
	cmd, err := NewData(1, alifn, treefn, cfreq)

	return cmd, err
}

// maxInt returns maximum integer value.
func maxInt(a int, b ...int) int {
	for _, v := range b {
		if v > a {
			a = v
		}
	}
	return a
}

// setLogLevel sets the default log-level to WARNING.
func setLogLevel() {
}

func strFltSlice(fs []float64) string {
	var b bytes.Buffer
	for _, f := range fs {
		fmt.Fprintf(&b, "%.3f ", f)
	}
	return b.String()
}

func floatRange(start, step float64, n uint) (res []float64) {
	res = make([]float64, n)
	for i := 0; i < int(n); i++ {
		res[i] = start
		start += step
	}
	return
}

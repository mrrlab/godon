package cmodel

import (
	"bytes"
	"fmt"
	"path"

	"github.com/op/go-logging"
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

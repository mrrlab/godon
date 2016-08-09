package cmodel

import (
	"os"
	"path"

	"github.com/op/go-logging"

	"bitbucket.org/Davydov/godon/bio"
	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/tree"
)

const (
	data1 = "EMGT00050000008747.Drosophila.002"
	data2 = "ENSGT00550000073950.Euteleostomi.07.001"
	data3 = "EMGT00050000000025.Drosophila.001"

	smallDiff = 1e-3
)

// log is a global logging variable.
var log = logging.MustGetLogger("cmodel")

// GetTreeAlignment returns a tree and elignment for testing purposes.
func GetTreeAlignment(data string) (t *tree.Tree, cali codon.CodonSequences, err error) {
	tf, err := os.Open(path.Join("testdata", data+".nwk"))
	if err != nil {
		return
	}
	defer tf.Close()

	t, err = tree.ParseNewick(tf)
	if err != nil {
		return
	}

	af, err := os.Open(path.Join("testdata", data+".fst"))
	if err != nil {
		return
	}
	defer af.Close()

	ali, err := bio.ParseFasta(af)
	if err != nil {
		return
	}

	cali, err = codon.ToCodonSequences(ali)
	if err != nil {
		return
	}

	return
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
	logging.SetLevel(logging.WARNING, "godon")
	logging.SetLevel(logging.WARNING, "optimize")
	logging.SetLevel(logging.WARNING, "cmodel")
}

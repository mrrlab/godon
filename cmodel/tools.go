package cmodel

import (
	"os"
	"path"

	"bitbucket.org/Davydov/godon/bio"
	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/tree"

	"github.com/op/go-logging"
)

const (
	data1 = "EMGT00050000008747.Drosophila.002"
	data2 = "ENSGT00550000073950.Euteleostomi.07.001"
	data3 = "EMGT00050000000025.Drosophila.001"

	smallDiff = 1e-4
)

var log = logging.MustGetLogger("cmodel")

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

func setLogLevel() {
	logging.SetLevel(logging.WARNING, "godon")
	logging.SetLevel(logging.WARNING, "optimize")
	logging.SetLevel(logging.WARNING, "cmodel")
}

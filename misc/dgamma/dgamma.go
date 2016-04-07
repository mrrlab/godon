package main

import (
	"flag"
	"fmt"

	"bitbucket.org/Davydov/godon/paml"
)

func main() {
	alpha := flag.Float64("alpha", 1, "alpha")
	ncat := flag.Int("ncat", 4, "ncat")
	useMedian := flag.Bool("median", false, "Use median instead of mean")
	flag.Parse()

	r := paml.DiscreteGamma(*alpha, *alpha, *ncat, *useMedian, nil, nil)
	fmt.Println(r)
}

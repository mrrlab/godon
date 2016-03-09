package main

import (
	"flag"
	"fmt"

	"bitbucket.org/Davydov/godon/paml"
)

func main() {
	p := flag.Float64("p", 1, "p")
	q := flag.Float64("q", 1, "q")
	ncat := flag.Int("ncat", 4, "ncat")
	useMedian := flag.Bool("median", false, "Use median instead of mean")
	flag.Parse()

	r := paml.DiscreteGamma(*p, *q, *ncat, *useMedian, nil, nil)
	fmt.Println(r)
}

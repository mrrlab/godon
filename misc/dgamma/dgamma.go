// DBeta returns values of a discre beta distribution.
package main

import (
	"flag"
	"fmt"

	"bitbucket.org/Davydov/godon/dist"
)

func main() {
	alpha := flag.Float64("alpha", 1, "alpha")
	ncat := flag.Int("ncat", 4, "ncat")
	useMedian := flag.Bool("median", false, "Use median instead of mean")
	flag.Parse()

	r := dist.DiscreteGamma(*alpha, *alpha, *ncat, *useMedian, nil, nil)
	fmt.Println(r)
}

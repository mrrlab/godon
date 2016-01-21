package main

import (
	"flag"
	"fmt"

	"github.com/gonum/plot"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/plotutil"
	"github.com/gonum/plot/vg"

	"bitbucket.com/Davydov/godon/paml"
)

func main() {
	alpha := flag.Float64("alpha", 1, "alpha")
	k := flag.Int("k", 4, "k")
	useMedian := flag.Bool("median", false, "Use median instead of mean")
	flag.Parse()

	a, b := paml.DiscreteGamma(*alpha, *alpha, *k, *useMedian)
	fmt.Println(a)
	fmt.Println(b)
	p, err := plot.New()
	if err != nil {
		panic(err)
	}

	pts := make(plotter.XYs, *k)
	x := 0.0
	for i, v := range b {
		pts[i].X = v
		pts[i].Y = x
		x += 1. / float64(*k)
	}

	err = plotutil.AddLinePoints(p,
		"first", pts)
	if err != nil {
		panic(err)
	}

	if err := p.Save(4*vg.Inch, 4*vg.Inch, "points.png"); err != nil {
		panic(err)
	}

}

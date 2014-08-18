package main

import (
	"math"
	"testing"
)

const SmallDiff = 1e-6

func TestReflect(tst *testing.T) {
	if math.Abs(Reflect(2, 0, 5)-2.0) > SmallDiff ||
		math.Abs(Reflect(2, -1, 5)-2.0) > SmallDiff ||
		math.Abs(Reflect(-2, 1, 5)-4.0) > SmallDiff ||
		math.Abs(Reflect(6, 1, 5)-4.0) > SmallDiff {
		tst.Error("Error reflecting")
	}
}

// Package dist implements functions for discrete distributions.
package dist

/*
this is largely based on PAML, which is licensed under GNU GPL v3.
Many functions where replaced by math/mathext equivalents,
DiscreteBeta now checks boundaries and behaves better.
*/

import (
	"math"

	"github.com/gonum/mathext"
)

/*

QuantileChi2 returns z so that Prob{x<z}=prob where x is Chi2
distributed with df=v

returns -1 if in error.  0.000002<prob<0.999998

RATNEST FORTRAN by Best DJ & Roberts DE (1975) The percentage points
of the Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)

Converted into C by Ziheng Yang, Oct. 1993.

*/
func QuantileChi2(prob, v float64) (ch float64) {
	e := .5e-6
	aa := .6931471805
	p := prob
	small := 1e-6
	a := 0.0
	q := 0.0
	p1 := 0.0
	p2 := 0.0
	t := 0.0
	x := 0.0
	b := 0.0
	//s1,s2,s3,s4,s5,s6;

	if p < small {
		return 0
	}
	if p > 1-small {
		return 9999
	}
	if v <= 0 {
		return -1
	}

	g, _ := math.Lgamma(v / 2)
	xx := v / 2
	c := xx - 1
	if v >= -1.24*math.Log(p) {
		goto l1
	}

	ch = math.Pow((p * xx * math.Exp(g+xx*aa)), 1/xx)
	if ch-e < 0 {
		return ch
	}
	goto l4
l1:
	if v > .32 {
		goto l3
	}
	ch = 0.4
	a = math.Log(1 - p)
l2:
	q = ch
	p1 = 1 + ch*(4.67+ch)
	p2 = ch * (6.73 + ch*(6.66+ch))
	t = -0.5 + (4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2
	ch -= (1 - math.Exp(a+g+.5*ch+c*aa)*p2/p1) / t
	if math.Abs(q/ch-1)-.01 <= 0 {
		goto l4
	} else {
		goto l2
	}
l3:
	x = QuantileNormal(p)
	p1 = 0.222222 / v
	ch = v * math.Pow((x*math.Sqrt(p1)+1-p1), 3.0)
	if ch > 2.2*v+6 {
		ch = -2 * (math.Log(1-p) - c*math.Log(.5*ch) + g)
	}
l4:
	q = ch
	p1 = .5 * ch
	t = IncompleteGamma(p1, xx)
	if t < 0 {
		panic("IncompleteGamma<0")
	}
	p2 = p - t
	t = p2 * math.Exp(xx*aa+g+p1-c*math.Log(ch))
	b = t / ch
	a = 0.5*t - b*c

	s1 := (210 + a*(140+a*(105+a*(84+a*(70+60*a))))) / 420
	s2 := (420 + a*(735+a*(966+a*(1141+1278*a)))) / 2520
	s3 := (210 + a*(462+a*(707+932*a))) / 2520
	s4 := (252 + a*(672+1182*a) + c*(294+a*(889+1740*a))) / 5040
	s5 := (84 + 264*a + c*(175+606*a)) / 2520
	s6 := (120 + c*(346+127*c)) / 5040
	ch += t * (1 + 0.5*t*s1 - b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))))
	if math.Abs(q/ch-1) > e {
		goto l4
	}

	return
}

// QuantileGamma returns quantile for gamma distribution.
func QuantileGamma(prob, alpha, beta float64) float64 {
	return QuantileChi2(prob, 2.0*(alpha)) / (2.0 * (beta))
}

// QuantileNormal returns quantile for normal distribution.
func QuantileNormal(prob float64) float64 {
	return mathext.NormalQuantile(prob)
}

/*

IncompleteGamma returns the incomplete gamma ratio I(x,alpha) where x
is the upper limit of the integration and alpha is the shape
parameter.

*/
func IncompleteGamma(x, alpha float64) (gin float64) {
	return mathext.GammaInc(alpha, x)
}

// DiscreteGamma returns discrete gamma distribution.
func DiscreteGamma(alpha, beta float64, K int, UseMedian bool, tmp, res []float64) []float64 {
	/*
	   discretization of G(alpha, beta) with equal proportions in each category.
	*/
	t := 0.0
	mean := alpha / beta

	if res == nil {
		res = make([]float64, K)
	}
	if tmp == nil {
		tmp = make([]float64, K)
	}

	if UseMedian { /* median */
		for i := 0; i < K; i++ {
			res[i] = QuantileGamma((float64(i)*2.+1)/(2.*float64(K)), alpha, beta)
		}
		for i := 0; i < K; i++ {
			t += res[i]
		}
		for i := 0; i < K; i++ {
			res[i] *= mean * float64(K) / t /* rescale so that the mean is alpha/beta. */
		}
	} else { /* mean */
		for i := 0; i < K-1; i++ { /* cutting points, Eq. 9 */
			tmp[i] = QuantileGamma((float64(i)+1.0)/float64(K), alpha, beta)
		}
		for i := 0; i < K-1; i++ { /* Eq. 10 */
			tmp[i] = IncompleteGamma(tmp[i]*beta, alpha+1)
		}
		res[0] = tmp[0] * mean * float64(K)
		for i := 1; i < K-1; i++ {
			res[i] = (tmp[i] - tmp[i-1]) * mean * float64(K)
		}
		res[K-1] = (1 - tmp[K-2]) * mean * float64(K)
	}

	return res
}

// LnBeta returns log of Beta function.
func LnBeta(p, q float64) float64 {
	lgp, _ := math.Lgamma(p)
	lgq, _ := math.Lgamma(q)
	lgpq, _ := math.Lgamma(p + q)
	return lgp + lgq - lgpq
}

var eps, alneps, sml, alnsml float64 = 0, 0, 0, 0

/*

CDFBeta returns distribution function of the standard form of the beta
distribution, that is, the incomplete beta ratio I_x(p,q).

This is also known as the incomplete beta function ratio I_x(p, q)

This is called from QuantileBeta() in a root-finding loop.

*/
func CDFBeta(x, pin, qin float64) float64 {
	return mathext.RegIncBeta(pin, qin, x)
}

/*
QuantileBeta calculates the Quantile of the beta distribution
*/
func QuantileBeta(prob, p, q float64) float64 {
	return mathext.InvRegIncBeta(p, q, prob)
}

// DiscreteBeta returns discrete beta distribution.
func DiscreteBeta(p, q float64, K int, UseMedian bool, tmp, res []float64) []float64 {
	/*
	   discretization of beta(p, q), with equal proportions in each category.
	*/
	mean := p / (p + q)
	t := 0.0

	if res == nil {
		res = make([]float64, K)
	}
	if tmp == nil {
		tmp = make([]float64, K)
	}

	if UseMedian { /* median */
		for i := 0; i < K; i++ {
			res[i] = QuantileBeta((float64(i)+0.5)/float64(K), p, q)
			t += res[i]
		}
		// normalization to keep the mean
		for i := 0; i < K; i++ {
			res[i] *= mean * float64(K) / t
		}
	} else { /* mean */
		for i := 0; i < K-1; i++ /* cutting points */ {
			tmp[i] = QuantileBeta((float64(i)+1.0)/float64(K), p, q)
		}
		tmp[K-1] = 1

		prevCdf := CDFBeta(tmp[0], p+1, q)

		res[0] = prevCdf * mean * float64(K)
		for i := 1; i < K; i++ { /* CDF */
			currCdf := CDFBeta(tmp[i], p+1, q)
			res[i] = (currCdf - prevCdf) * mean * float64(K)
			prevCdf = currCdf
		}

		for i := 0; i < K; i++ { /* correct out of region */
			lower := 0.0
			upper := tmp[i]
			if i > 0 {
				lower = tmp[i-1]
			}
			if res[i] < lower || res[i] > upper {
				// if upper-lower > 1e-3 {
				// 	fmt.Println("out of bounds", res[i], lower, upper)
				// 	fmt.Println(p, q, K, UseMedian)
				// }

				// switch to median
				res[i] = QuantileBeta((float64(i)+0.5)/float64(K), p, q)
				if res[i] < lower || res[i] > upper { //out of bounds again
					// if upper-lower > 1e-3 {
					// 	fmt.Println("out of bounds, again")
					// }

					// switch to average
					res[i] = (upper + lower) / 2
				}
			}
		}
	}

	return res
}

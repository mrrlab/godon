package paml

/* this code comes from PAML 4.8a */

import "math"

func QuantileChi2(prob, v float64) (ch float64) {
	/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
	   returns -1 if in error.   0.000002<prob<0.999998
	   RATNEST FORTRAN by
	       Best DJ & Roberts DE (1975) The percentage points of the
	       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
	   Converted into C by Ziheng Yang, Oct. 1993.
	*/
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
	t = IncompleteGamma(p1, xx, g)
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

func QuantileGamma(prob, alpha, beta float64) float64 {
	return QuantileChi2(prob, 2.0*(alpha)) / (2.0 * (beta))
}

func QuantileNormal(prob float64) float64 {
	/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
	   returns (-9999) if in error
	   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
	   Applied Statistics 22: 96-97 (AS70)

	   Newer methods:
	     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
	       normal distribution.  37: 477-484.
	     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
	       points of the normal distribution.  26: 118-121.
	*/
	a0 := -.322232431088
	a1 := -1.0
	a2 := -.342242088547
	a3 := -.0204231210245
	a4 := -.453642210148e-4
	b0 := .0993484626060
	b1 := .588581570495
	b2 := .531103462366
	b3 := .103537752850
	b4 := .0038560700634
	var z, p1 float64
	p := prob

	if p < 0.5 {
		p1 = p
	} else {
		p1 = 1 - p
	}
	if p1 < 1e-20 {
		z = 999
	} else {
		y := math.Sqrt(math.Log(1 / (p1 * p1)))
		z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0)/((((y*b4+b3)*y+b2)*y+b1)*y+b0)
	}
	if p < 0.5 {
		return -z
	} else {
		return z
	}
}

func IncompleteGamma(x, alpha, ln_gamma_alpha float64) (gin float64) {
	/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
	           limit of the integration and alpha is the shape parameter.
	   returns (-1) if in error
	   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
	   (1) series expansion,     if (alpha>x || x<=1)
	   (2) continued fraction,   otherwise
	   RATNEST FORTRAN by
	   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
	   19: 285-287 (AS32)
	*/
	i := 0
	p := alpha
	g := ln_gamma_alpha
	accurate := 1e-10
	overflow := 1e60
	rn := 0.0
	a := 0.0
	b := 0.0
	an := 0.0
	dif := 0.0
	term := 0.0
	pn := make([]float64, 6)

	if x == 0 {
		return 0
	}
	if x < 0 || p <= 0 {
		return -1
	}

	factor := math.Exp(p*math.Log(x) - x - g)
	if x > 1 && x >= p {
		goto l30
	}
	/* (1) series expansion */
	gin = 1.0
	term = 1.0
	rn = p
l20:
	rn++
	term *= x / rn
	gin += term
	if term > accurate {
		goto l20
	}
	gin *= factor / p
	goto l50
l30:
	/* (2) continued fraction */
	a = 1 - p
	b = a + x + 1
	term = 0
	pn[0] = 1
	pn[1] = x
	pn[2] = x + 1
	pn[3] = x * b
	gin = pn[2] / pn[3]
l32:
	a++
	b += 2
	term++
	an = a * term
	for i = 0; i < 2; i++ {
		pn[i+4] = b*pn[i+2] - an*pn[i]
	}
	if pn[5] == 0 {
		goto l35
	}
	rn = pn[4] / pn[5]
	dif = math.Abs(gin - rn)
	if dif > accurate {
		goto l34
	}
	if dif <= accurate*rn {
		goto l42
	}
l34:
	gin = rn
l35:
	for i = 0; i < 4; i++ {
		pn[i] = pn[i+2]
	}
	if math.Abs(pn[4]) < overflow {
		goto l32
	}
	for i = 0; i < 4; i++ {
		pn[i] /= overflow
	}
	goto l32
l42:
	gin = 1 - factor*gin

l50:
	return
}

func DiscreteGamma(alpha, beta float64, K int, UseMedian bool, freqK, rK []float64) []float64 {
	/*
	   discretization of G(alpha, beta) with equal proportions in each category.
	*/
	i := 0
	t := 0.0
	mean := alpha / beta

	if rK == nil {
		rK = make([]float64, K)
	}
	if freqK == nil {
		freqK = make([]float64, K)
	}

	if UseMedian { /* median */
		for i = 0; i < K; i++ {
			rK[i] = QuantileGamma((float64(i)*2.+1)/(2.*float64(K)), alpha, beta)
		}
		for i = 0; i < K; i++ {
			t += rK[i]
		}
		for i = 0; i < K; i++ {
			rK[i] *= mean * float64(K) / t /* rescale so that the mean is alpha/beta. */
		}
	} else { /* mean */
		lnga1, _ := math.Lgamma(alpha + 1)
		for i = 0; i < K-1; i++ { /* cutting points, Eq. 9 */
			freqK[i] = QuantileGamma((float64(i)+1.0)/float64(K), alpha, beta)
		}
		for i = 0; i < K-1; i++ { /* Eq. 10 */
			freqK[i] = IncompleteGamma(freqK[i]*beta, alpha+1, lnga1)
		}
		rK[0] = freqK[0] * mean * float64(K)
		for i = 1; i < K-1; i++ {
			rK[i] = (freqK[i] - freqK[i-1]) * mean * float64(K)
		}
		rK[K-1] = (1 - freqK[K-2]) * mean * float64(K)
	}

	return rK
}

package paml

/* this code comes from PAML 4.8a */

import (
	"fmt"
	"math"
)

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
		lnga1, _ := math.Lgamma(alpha + 1)
		for i := 0; i < K-1; i++ { /* cutting points, Eq. 9 */
			tmp[i] = QuantileGamma((float64(i)+1.0)/float64(K), alpha, beta)
		}
		for i := 0; i < K-1; i++ { /* Eq. 10 */
			tmp[i] = IncompleteGamma(tmp[i]*beta, alpha+1, lnga1)
		}
		res[0] = tmp[0] * mean * float64(K)
		for i := 1; i < K-1; i++ {
			res[i] = (tmp[i] - tmp[i-1]) * mean * float64(K)
		}
		res[K-1] = (1 - tmp[K-2]) * mean * float64(K)
	}

	return res
}

func LnBeta(p, q float64) float64 {
	lgp, _ := math.Lgamma(p)
	lgq, _ := math.Lgamma(q)
	lgpq, _ := math.Lgamma(p + q)
	return lgp + lgq - lgpq
}

var eps, alneps, sml, alnsml float64 = 0, 0, 0, 0

func CDFBeta(x, pin, qin, lnbeta float64) float64 {
	/* Returns distribution function of the standard form of the beta distribution,
	   that is, the incomplete beta ratio I_x(p,q).

	   This is also known as the incomplete beta function ratio I_x(p, q)

	   lnbeta is log of the complete beta function; provide it if known,
	   and otherwise use 0.

	   This is called from QuantileBeta() in a root-finding loop.

	    This routine is a translation into C of a Fortran subroutine
	    by W. Fullerton of Los Alamos Scientific Laboratory.
	    Bosten and Battiste (1974).
	    Remark on Algorithm 179, CACM 17, p153, (1974).
	*/
	var ans, c, finsum, p, ps, p1, q, term, xb, xi, y float64
	small := 1e-15
	var n, ib int

	if x < small {
		return 0
	}
	if x > 1-small {
		return 1
	}

	if pin <= 0 || qin <= 0 {
		panic(fmt.Sprintf("p=%.4f q=%.4f: parameter outside range in CDFBeta", pin, qin))
	}

	if eps == 0 { /* initialize machine constants ONCE */
		eps = 1.203e-16
		alneps = math.Log(eps)
		sml = 2.22508e-308
		alnsml = math.Log(sml)
	}
	y = x
	p = pin
	q = qin

	/* swap tails if x is greater than the mean */
	if p/(p+q) < x {
		y = 1 - y
		p = qin
		q = pin
	}

	if lnbeta == 0 {
		lnbeta = LnBeta(p, q)
	}

	if (p+q)*y/(p+1) < eps { /* tail approximation */
		ans = 0
		xb = p*math.Log(math.Max(y, sml)) - math.Log(p) - lnbeta
		if xb > alnsml && y != 0 {
			ans = math.Exp(xb)
		}
		if y != x || p != pin {
			ans = 1 - ans
		}
	} else {
		/* evaluate the infinite sum first.  term will equal */
		/* y^p / beta(ps, p) * (1 - ps)-sub-i * y^i / fac(i) */
		ps = q - math.Floor(q)
		if ps == 0 {
			ps = 1
		}

		xb = LnBeta(ps, p)
		xb = p*math.Log(y) - xb - math.Log(p)

		ans = 0
		if xb >= alnsml {
			ans = math.Exp(xb)
			term = ans * p
			if ps != 1 {
				n = int(math.Max(alneps/math.Log(y), 4.0))
			}
			for i := 1; i <= n; i++ {
				xi = float64(i)
				term = term * (xi - ps) * y / xi
				ans = ans + term/(p+xi)
			}
		}

		/* evaluate the finite sum. */
		if q > 1 {
			xb = p*math.Log(y) + q*math.Log(1-y) - lnbeta - math.Log(q)
			ib = (int)(xb / alnsml)
			if ib < 0 {
				ib = 0
			}
			term = math.Exp(xb - float64(ib)*alnsml)
			c = 1 / (1 - y)
			p1 = q * c / (p + q - 1)

			finsum = 0
			n = int(q)
			if q == float64(n) {
				n = n - 1
			}
			for i := 1; i <= n; i++ {
				if p1 <= 1 && term/eps <= finsum {
					break
				}
				xi = float64(i)
				term = (q - xi + 1) * c * term / (p + q - xi)
				if term > 1 {
					ib = ib - 1
					term = term * sml
				}
				if ib == 0 {
					finsum = finsum + term
				}
			}
			ans = ans + finsum
		}
		if y != x || p != pin {
			ans = 1 - ans
		}
		if ans > 1 {
			ans = 1
		}
		if ans < 0 {
			ans = 0
		}
	}
	return ans
}

func QuantileBeta(prob, p, q, lnbeta float64) float64 {
	/* This calculates the Quantile of the beta distribution

	   Cran, G. W., K. J. Martin and G. E. Thomas (1977).
	   Remark AS R19 and Algorithm AS 109, Applied Statistics, 26(1), 111-114.
	   Remark AS R83 (v.39, 309-310) and correction (v.40(1) p.236).

	   My own implementation of the algorithm did not bracket the variable well.
	   This version is Adpated from the pbeta and qbeta routines from
	   "R : A Computer Language for Statistical Data Analysis".  It fails for
	   extreme values of p and q as well, although it seems better than my
	   previous version.
	   Ziheng Yang, May 2001
	*/
	fpu := 3e-308
	acu_min := 1e-300
	lower := fpu
	upper := 1 - 2.22e-16
	/* acu_min>= fpu: Minimal value for accuracy 'acu' which will depend on (a,p); */
	var i_pb, i_inn int
	var swap_tail bool
	niterations := 2000
	var a, adj, g, h, pp, prev, qq, r, s, t, tx, w, y, yprev float64
	var acu, xinbta float64

	if prob < 0 || prob > 1 || p < 0 || q < 0 {
		panic("out of range in QuantileBeta")
	}

	/* define accuracy and initialize */
	xinbta = prob

	/* test for admissibility of parameters */
	if p < 0 || q < 0 || prob < 0 || prob > 1 {
		panic("beta par err")
	}
	if prob == 0 || prob == 1 {
		return prob
	}

	if lnbeta == 0 {
		lnbeta = LnBeta(p, q)
	}

	/* change tail if necessary;  afterwards   0 < a <= 1/2    */
	if prob <= 0.5 {
		a = prob
		pp = p
		qq = q
		swap_tail = false
	} else {
		a = 1. - prob
		pp = q
		qq = p
		swap_tail = true
	}

	/* calculate the initial approximation */
	r = math.Sqrt(-math.Log(a * a))
	y = r - (2.30753+0.27061*r)/(1.+(0.99229+0.04481*r)*r)

	if pp > 1. && qq > 1. {
		r = (y*y - 3.) / 6.
		s = 1. / (pp*2. - 1.)
		t = 1. / (qq*2. - 1.)
		h = 2. / (s + t)
		w = y*math.Sqrt(h+r)/h - (t-s)*(r+5./6.-2./(3.*h))
		xinbta = pp / (pp + qq*math.Exp(w+w))
	} else {
		r = qq * 2.
		t = 1. / (9. * qq)
		t = r * math.Pow(1.-t+y*math.Sqrt(t), 3.)
		if t <= 0. {
			xinbta = 1. - math.Exp((math.Log((1.-a)*qq)+lnbeta)/qq)
		} else {
			t = (4.*pp + r - 2.) / t
			if t <= 1. {
				xinbta = math.Exp((math.Log(a*pp) + lnbeta) / pp)
			} else {
				xinbta = 1. - 2./(t+1.)
			}
		}
	}

	/* solve for x by a modified newton-raphson method, using CDFBeta */
	r = 1. - pp
	t = 1. - qq
	yprev = 0.
	adj = 1.

	/* Changes made by Ziheng to fix a bug in qbeta()
	   qbeta(0.25, 0.143891, 0.05) = 3e-308   wrong (correct value is 0.457227)
	*/
	if xinbta <= lower || xinbta >= upper {
		xinbta = (a + .5) / 2
	}

	/* Desired accuracy should depend on (a,p)
	 * This is from Remark .. on AS 109, adapted.
	 * However, it's not clear if this is "optimal" for IEEE double prec.
	 * acu = fmax2(acu_min, pow(10., -25. - 5./(pp * pp) - 1./(a * a)));
	 * NEW: 'acu' accuracy NOT for squared adjustment, but simple;
	 * ---- i.e.,  "new acu" = sqrt(old acu)
	 */
	acu = math.Pow(10., -13.-2.5/(pp*pp)-0.5/(a*a))
	acu = math.Max(acu, acu_min)

	for i_pb = 0; i_pb < niterations; i_pb++ {
		y = CDFBeta(xinbta, pp, qq, lnbeta)
		y = (y - a) * math.Exp(lnbeta+r*math.Log(xinbta)+t*math.Log(1.-xinbta))
		if y*yprev <= 0 {
			prev = math.Max(math.Abs(adj), fpu)
		}
		g = 1
		for i_inn = 0; i_inn < niterations; i_inn++ {
			adj = g * y
			if math.Abs(adj) < prev {
				tx = xinbta - adj /* trial new x */
				if tx >= 0. && tx <= 1. {
					if prev <= acu || math.Abs(y) <= acu {
						goto L_converged
					}
					if tx != 0. && tx != 1. {
						break
					}
				}
			}
			g /= 3.
		}
		if math.Abs(tx-xinbta) < fpu {
			goto L_converged
		}
		xinbta = tx
		yprev = y
	}

L_converged:
	if swap_tail {
		return 1. - xinbta
	}
	return xinbta
}

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

	lnbeta := LnBeta(p, q)
	if UseMedian { /* median */
		for i := 0; i < K; i++ {
			res[i] = QuantileBeta((float64(i)+0.5)/float64(K), p, q, lnbeta)
			t += res[i]
		}
		for i := 0; i < K; i++ {
			res[i] *= mean * float64(K) / t
		}
	} else { /* mean */
		for i := 0; i < K-1; i++ /* cutting points */ {
			tmp[i] = QuantileBeta((float64(i)+1.0)/float64(K), p, q, lnbeta)
		}
		tmp[K-1] = 1

		lnbeta1 := lnbeta - math.Log(1+q/p)
		for i := 0; i < K-1; i++ { /* CDF */
			tmp[i] = CDFBeta(tmp[i], p+1, q, lnbeta1)
		}
		res[0] = tmp[0] * mean * float64(K)
		for i := 1; i < K-1; i++ {
			res[i] = (tmp[i] - tmp[i-1]) * mean * float64(K)
		}
		res[K-1] = (1 - tmp[K-2]) * mean * float64(K)

		for i := 0; i < K; i++ {
			t += res[i] / float64(K)
		}
	}

	return res
}

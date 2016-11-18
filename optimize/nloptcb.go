package optimize

/*
#include <nlopt.h>
#include <math.h>
*/
import "C"

import (
	"math"
	"unsafe"
)

//export nloptCallback
func nloptCallback(n uint, x *C.double, grad *C.double, fData unsafe.Pointer) C.double {
	nlopt := lookupObject(uintptr(fData)).(*NLOPT)

	select {
	case s := <-nlopt.sig:
		log.Warningf("Received signal (%v), exiting.", s)
		C.nlopt_force_stop(nlopt.gopt)
		nlopt.stop = true
	default:
	}

	if nlopt.stop {
		return (C.double)(math.Inf(-1))
	}

	xsl := (*[1 << 30]C.double)(unsafe.Pointer(x))[:n:n]

	for i := range xsl {
		nlopt.parameters[i].Set((float64)(xsl[i]))
	}
	l1 := nlopt.Likelihood()
	nlopt.calls++
	if grad != nil {
		gradsl := (*[1 << 30]C.double)(unsafe.Pointer(grad))[:n:n]
	g:
		for i := range xsl {
			inv := false
			v := (float64)(xsl[i]) + nlopt.dH

			// this shouldn't happen with current boundaries
			// but to be safe
			if v >= nlopt.parameters[i].GetMax() {
				v = (float64)(xsl[i]) - nlopt.dH
				inv = true
			}

			nlopt.parameters[i].Set(v)
			l2 := nlopt.Likelihood()
			nlopt.calls++

			gradsl[i] = (C.double)((l2 - l1) / nlopt.dH)
			if inv {
				gradsl[i] = -gradsl[i]
			}

			nlopt.parameters[i].Set((float64)(xsl[i]))

			select {
			case s := <-nlopt.sig:
				log.Warningf("Received signal (%v), exiting.", s)
				C.nlopt_force_stop(nlopt.gopt)
				nlopt.stop = true
				break g
			default:
			}
		}
	}

	if !nlopt.stop {
		if nlopt.i%nlopt.repPeriod == 0 {
			nlopt.PrintLine(nlopt.parameters, l1)
		}
		nlopt.i++
	}

	return (C.double)(l1)
}

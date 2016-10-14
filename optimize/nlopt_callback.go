package optimize

/*
#include <nlopt.h>
#include <math.h>
*/
import "C"

import (
	"math"
	"sync"
	"unsafe"
)

// objects is a container for the actual go objects. We cannot pass a
// raw pointer to go object if it has pointer inside, so we can pass
// its' index instead.
var objects = make(map[int]interface{})

// objectIndex stores an index to use for new object.
var objectIndex int

// objectMutex is a mutex preventing simultanious access to the object
// storage.
var objectMutex sync.Mutex

// registerObject registers a new object and returns its' index
// (>=1).
func registerObject(obj interface{}) int {
	objectMutex.Lock()
	defer objectMutex.Unlock()
	// ensure minimum index is 1.
	objectIndex++
	for objects[objectIndex] != nil {
		objectIndex++
	}
	objects[objectIndex] = obj
	return objectIndex
}

// lookupCallback returns an object given an index.
func lookupObject(i int) interface{} {
	objectMutex.Lock()
	defer objectMutex.Unlock()
	return objects[i]
}

// unregisterObject unregisters an object  by removing it from the
// objects map.
func unregisterObject(i int) {
	objectMutex.Lock()
	defer objectMutex.Unlock()
	delete(objects, i)
}

//export nlopt_callback
func nlopt_callback(n uint, x *C.double, grad *C.double, f_data unsafe.Pointer) C.double {
	nlopt := lookupObject(*(*int)(f_data)).(*NLOPT)

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

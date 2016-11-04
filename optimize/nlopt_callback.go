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
var objects = make(map[uint]interface{})

// objectIndex stores an index to use for new object.
var objectIndex uint

// objectMutex is a mutex preventing simultanious access to the object
// storage.
var objectMutex sync.Mutex

// registerObject registers a new object and returns its' index
// (>=1).
func registerObject(obj interface{}) uint {
	objectMutex.Lock()
	defer objectMutex.Unlock()
	// We always increment objectIndex to have more or less unique
	// ids. This way it is easier to debug problems with reusing
	// unregistered ids.
	objectIndex++
	startIndex := objectIndex
	for objectIndex == 0 || objects[objectIndex] != nil {
		objectIndex++
		// If the map is full, i.e. all non-zero uints were
		// used, we do not want to loop infinitely. We check
		// if we already encountered the starting index. If
		// so, we panic. In practice this is very unlikely to
		// have this kind of problem since all the objects are
		// unregistered at the end of the function call.
		if objectIndex == startIndex {
			panic("no more space in the map to store an object")
		}

	}
	log.Debugf("register: %v (%p)", objectIndex, obj)
	objects[objectIndex] = obj
	return objectIndex
}

// lookupCallback returns an object given an index.
func lookupObject(i uint) interface{} {
	objectMutex.Lock()
	defer objectMutex.Unlock()
	return objects[i]
}

// unregisterObject unregisters an object  by removing it from the
// objects map.
func unregisterObject(i uint) {
	objectMutex.Lock()
	log.Debugf("unregister: %v (%p)", i, objects[i])
	defer objectMutex.Unlock()
	delete(objects, i)
}

//export nloptCallback
func nloptCallback(n uint, x *C.double, grad *C.double, fData unsafe.Pointer) C.double {
	nlopt := lookupObject(*(*uint)(fData)).(*NLOPT)

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

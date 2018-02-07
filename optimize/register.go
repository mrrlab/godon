package optimize

import "sync"

// objects is a container for the actual go objects. We cannot pass a
// raw pointer to go object if it has pointer inside, so we can pass
// its' index instead.
var objects = make(map[uintptr]interface{})

// objectIndex stores an index to use for new object.
var objectIndex uintptr

// objectMutex is a mutex preventing simultaneous access to the object
// storage.
var objectMutex sync.Mutex

// registerObject registers a new object and returns its' index
// (>=1).
func registerObject(obj interface{}) uintptr {
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
func lookupObject(i uintptr) interface{} {
	objectMutex.Lock()
	defer objectMutex.Unlock()
	return objects[i]
}

// unregisterObject unregisters an object  by removing it from the
// objects map.
func unregisterObject(i uintptr) {
	objectMutex.Lock()
	log.Debugf("unregister: %v (%p)", i, objects[i])
	defer objectMutex.Unlock()
	delete(objects, i)
}

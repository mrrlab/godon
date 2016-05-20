package optimize

// #cgo LDFLAGS: -lnlopt -lm
/*
#include <nlopt.h>

extern double nlopt_callback(unsigned int, const double *x, double *grad, void *f_data);


double callback_adaptor(unsigned n, const double *x, double *grad, void *f_data) {
  return nlopt_callback(n, x, grad, f_data);
}
*/
import "C"

import (
	"unsafe"
)

const (
	NLOPT_COBYLA = iota
	NLOPT_BOBYQA
	NLOPT_SIMPLEX
	NLOPT_LBFGS
	NLOPT_SQP
)

type NLOPT struct {
	BaseOptimizer
	gopt C.nlopt_opt
	//lopt C.nlopt_opt
	dH          float64
	stop        bool
	finishLBFGS bool
	seed        int64
	xtol_rel    float64
	ftol_rel    float64
	algorithm   C.nlopt_algorithm
}

func NewNLOPT(algorithm int, seed int64) (nlopt *NLOPT) {
	nlopt = &NLOPT{
		BaseOptimizer: BaseOptimizer{
			repPeriod: 10,
		},
		dH:          1e-6,
		finishLBFGS: true,
		seed:        seed,
		xtol_rel:    1e-9,
		ftol_rel:    1e-9,
	}
	switch algorithm {
	case NLOPT_COBYLA:
		nlopt.algorithm = C.NLOPT_LN_COBYLA
	case NLOPT_BOBYQA:
		nlopt.algorithm = C.NLOPT_LN_BOBYQA
	case NLOPT_SIMPLEX:
		nlopt.algorithm = C.NLOPT_LN_NELDERMEAD
	case NLOPT_LBFGS:
		nlopt.algorithm = C.NLOPT_LD_LBFGS
	case NLOPT_SQP:
		nlopt.algorithm = C.NLOPT_LD_SLSQP
	default:
		log.Fatalf("Unknown algorithm specified (%d).", algorithm)
	}
	log.Infof("NLopt algorithm: %v.", C.GoString(C.nlopt_algorithm_name(nlopt.algorithm)))
	return
}

func (n *NLOPT) Run(iterations int) {

	n.gopt = C.nlopt_create(n.algorithm, (C.uint)(len(n.parameters)))
	defer C.nlopt_destroy(n.gopt)
	//n.lopt = C.nlopt_create(C.NLOPT_LD_LBFGS, (C.uint)(len(n.parameters)))
	//defer C.nlopt_destroy(n.lopt)
	//C.nlopt_set_ftol_rel(n.lopt, 1e-3)
	//C.nlopt_set_xtol_rel(n.lopt, 1e-5)
	//C.nlopt_set_local_optimizer(n.gopt, n.lopt)

	C.nlopt_set_max_objective(n.gopt, (C.nlopt_func)(unsafe.Pointer(C.callback_adaptor)), unsafe.Pointer(n))

	lb := make([]C.double, len(n.parameters))
	ub := make([]C.double, len(n.parameters))
	x := make([]C.double, len(n.parameters))
	for i, par := range n.parameters {
		lb[i] = (C.double)(par.GetMin())
		ub[i] = (C.double)(par.GetMax())
		x[i] = (C.double)(par.Get())
	}
	C.nlopt_set_lower_bounds(n.gopt, &lb[0])
	C.nlopt_set_upper_bounds(n.gopt, &ub[0])

	C.nlopt_set_xtol_rel(n.gopt, (C.double)(n.xtol_rel))
	C.nlopt_set_ftol_rel(n.gopt, (C.double)(n.ftol_rel))

	var maxf C.double

	if n.seed > 0 {
		C.nlopt_srand((C.ulong)(n.seed))
	}

	n.PrintHeader()
	res := C.nlopt_optimize(n.gopt, (*C.double)(unsafe.Pointer(&x[0])), &maxf)
	if res < 0 {
		log.Fatalf("nlopt failed with code: %v", res)
	}

	for i, par := range n.parameters {
		par.Set((float64)(x[i]))
	}

	n.maxL = (float64)(maxf)
	n.maxLPar = n.parameters.Values(n.maxLPar)
}

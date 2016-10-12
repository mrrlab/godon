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

// NLopt optimization methods.
const (
	// Constrained Optimization BY Linear Approximations
	// (local, gradient-free)
	NLOPT_COBYLA = iota
	// BOBYQA optimizer (local, gradient-free)
	NLOPT_BOBYQA
	// Downhill simplex (local, gradient-free)
	NLOPT_SIMPLEX
	// LBFGS: local Broyden–Fletcher–Goldfarb–Shanno
	// (local, gradient-based)
	NLOPT_LBFGS
	// SQP: sequential quadratic programming
	// (local, gradient-based)
	NLOPT_SQP
	// DIviding RECTangles algorithm (global)
	NLOPT_DIRECT
	// Controlled Random Search (global)
	NLOPT_CRS
	// Multi-Level Single-Linkage (global)
	NLOPT_MLSL
)

const (
	// How much to shrink [min; max] interval.
	shrink = 1e-7
)

// returnStatus converts status to a constant name.
var returnStatus = map[C.nlopt_result]string{
	1:  "NLOPT_SUCCESS",
	2:  "NLOPT_STOPVAL_REACHED",
	3:  "NLOPT_FTOL_REACHED",
	4:  "NLOPT_XTOL_REACHED",
	5:  "NLOPT_MAXEVAL_REACHED",
	6:  "NLOPT_MAXTIME_REACHED",
	-1: "NLOPT_FAILURE",
	-2: "NLOPT_INVALID_ARGS",
	-3: "NLOPT_OUT_OF_MEMORY",
	-4: "NLOPT_ROUNDOFF_LIMITED",
	-5: "NLOPT_FORCED_STOP",
}

// NLOPT is nlopt optimizer.
type NLOPT struct {
	BaseOptimizer
	gopt          C.nlopt_opt
	lopt          C.nlopt_opt
	dH            float64
	stop          bool
	finishLBFGS   bool
	seed          int64
	ftol_rel      float64
	ftol_abs      float64
	xtol_rel      float64
	algorithm     C.nlopt_algorithm
	loc_ftol_rel  float64
	loc_ftol_abs  float64
	loc_xtol_rel  float64
	loc_algorithm C.nlopt_algorithm
}

// NewNLOPT creates a new NLOPT optimizer.
func NewNLOPT(algorithm int, seed int64) (nlopt *NLOPT) {
	nlopt = &NLOPT{
		BaseOptimizer: BaseOptimizer{
			repPeriod: 10,
		},
		dH:           1e-6,
		finishLBFGS:  true,
		seed:         seed,
		ftol_rel:     1e-9,
		ftol_abs:     1e-4,
		xtol_rel:     1e-5,
		loc_ftol_rel: 1e-8,
		loc_ftol_abs: 1e-4,
		loc_xtol_rel: 1e-2,
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
	case NLOPT_DIRECT:
		nlopt.algorithm = C.NLOPT_GN_DIRECT_L_RAND
		nlopt.ftol_rel = 1e-11
		nlopt.ftol_abs = 1e-7
		nlopt.xtol_rel = 1e-4
	case NLOPT_CRS:
		nlopt.algorithm = C.NLOPT_GN_CRS2_LM
		nlopt.ftol_rel = 1e-7
		nlopt.ftol_abs = 1e-3
		nlopt.xtol_rel = 1e-4
	case NLOPT_MLSL:
		nlopt.algorithm = C.NLOPT_GN_MLSL_LDS
		nlopt.ftol_rel = 1e-1
		nlopt.ftol_abs = 1e-1
		nlopt.xtol_rel = 1e-1
		nlopt.loc_algorithm = C.NLOPT_LN_BOBYQA
		nlopt.loc_ftol_rel = 1e-4
		nlopt.loc_ftol_abs = 1
		nlopt.loc_xtol_rel = 1e-2
	default:
		log.Fatalf("Unknown algorithm specified (%d).", algorithm)
	}
	return
}

// Run runs the optimizer.
func (n *NLOPT) Run(iterations int) {
	log.Infof("NLopt algorithm: %v.", C.GoString(C.nlopt_algorithm_name(n.algorithm)))
	n.gopt = C.nlopt_create(n.algorithm, (C.uint)(len(n.parameters)))
	defer C.nlopt_destroy(n.gopt)
	if n.loc_algorithm != 0 {
		log.Infof("local algorithm: %v.", C.GoString(C.nlopt_algorithm_name(n.loc_algorithm)))
		n.lopt = C.nlopt_create(n.loc_algorithm, (C.uint)(len(n.parameters)))
		defer C.nlopt_destroy(n.lopt)
		C.nlopt_set_population(n.gopt, 1)
		C.nlopt_set_ftol_rel(n.lopt, (C.double)(n.loc_ftol_rel))
		C.nlopt_set_ftol_abs(n.lopt, (C.double)(n.loc_ftol_abs))
		C.nlopt_set_xtol_rel(n.lopt, (C.double)(n.loc_xtol_rel))
		C.nlopt_set_local_optimizer(n.gopt, n.lopt)
	}

	C.nlopt_set_max_objective(n.gopt, (C.nlopt_func)(unsafe.Pointer(C.callback_adaptor)), unsafe.Pointer(n))

	lb := make([]C.double, len(n.parameters))
	ub := make([]C.double, len(n.parameters))
	x := make([]C.double, len(n.parameters))
	for i, par := range n.parameters {
		val := par.Get()
		nMin := par.GetMin() + shrink
		nMax := par.GetMax() - shrink
		if val > nMax {
			val = nMax
		}
		if val < nMin {
			val = nMin
		}
		par.Set(val)
		lb[i] = (C.double)(nMin)
		ub[i] = (C.double)(nMax)
		x[i] = (C.double)(val)
	}
	C.nlopt_set_lower_bounds(n.gopt, &lb[0])
	C.nlopt_set_upper_bounds(n.gopt, &ub[0])

	log.Infof("Stopping criteria: ftol_rel=%g, ftol_abs=%g, xtol_rel=%g, maxeval=%d",
		n.ftol_rel, n.ftol_abs, n.xtol_rel, iterations)
	C.nlopt_set_ftol_rel(n.gopt, (C.double)(n.ftol_rel))
	C.nlopt_set_ftol_abs(n.gopt, (C.double)(n.ftol_abs))
	C.nlopt_set_xtol_rel(n.gopt, (C.double)(n.xtol_rel))
	C.nlopt_set_maxeval(n.gopt, (C.int)(iterations))

	var maxf C.double

	if n.seed > 0 {
		C.nlopt_srand((C.ulong)(n.seed))
	}

	n.PrintHeader()
	res := C.nlopt_optimize(n.gopt, (*C.double)(unsafe.Pointer(&x[0])), &maxf)
	if res < 0 {
		log.Fatalf("nlopt failed with code: %v (%v)", res, returnStatus[res])
	} else {
		log.Infof("nlopt success with code: %d (%v)", res, returnStatus[res])
	}

	for i, par := range n.parameters {
		par.Set((float64)(x[i]))
	}

	n.maxL = (float64)(maxf)
	n.maxLPar = n.parameters.Values(n.maxLPar)
}

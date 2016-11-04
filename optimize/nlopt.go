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

// returnStatus converts status to a constant name.
var returnStatus = map[C.nlopt_result]string{
	// 0 is the default value, if there's no result yet
	0:  "NO_RESULT",
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
	gopt         C.nlopt_opt
	lopt         C.nlopt_opt
	dH           float64
	stop         bool
	finishLBFGS  bool
	seed         int64
	ftolRel      float64
	ftolAbs      float64
	xtolRel      float64
	algorithm    C.nlopt_algorithm
	locFtolRel   float64
	locFtolAbs   float64
	locXtolRel   float64
	optRes       C.nlopt_result
	locAlgorithm C.nlopt_algorithm
}

// NewNLOPT creates a new NLOPT optimizer.
func NewNLOPT(algorithm int, seed int64) (nlopt *NLOPT) {
	nlopt = &NLOPT{
		BaseOptimizer: BaseOptimizer{
			repPeriod: 10,
		},
		dH:          1e-6,
		finishLBFGS: true,
		seed:        seed,
		ftolRel:     1e-9,
		ftolAbs:     1e-4,
		xtolRel:     1e-5,
		locFtolRel:  1e-8,
		locFtolAbs:  1e-4,
		locXtolRel:  1e-2,
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
		nlopt.ftolRel = 1e-11
		nlopt.ftolAbs = 1e-7
		nlopt.xtolRel = 1e-4
	case NLOPT_CRS:
		nlopt.algorithm = C.NLOPT_GN_CRS2_LM
		nlopt.ftolRel = 1e-7
		nlopt.ftolAbs = 1e-3
		nlopt.xtolRel = 1e-4
	case NLOPT_MLSL:
		nlopt.algorithm = C.NLOPT_GN_MLSL_LDS
		nlopt.ftolRel = 1e-1
		nlopt.ftolAbs = 1e-1
		nlopt.xtolRel = 1e-1
		nlopt.locAlgorithm = C.NLOPT_LN_BOBYQA
		nlopt.locFtolRel = 1e-4
		nlopt.locFtolAbs = 1
		nlopt.locXtolRel = 1e-2
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
	if n.locAlgorithm != 0 {
		log.Infof("local algorithm: %v.", C.GoString(C.nlopt_algorithm_name(n.locAlgorithm)))
		n.lopt = C.nlopt_create(n.locAlgorithm, (C.uint)(len(n.parameters)))
		defer C.nlopt_destroy(n.lopt)
		C.nlopt_set_population(n.gopt, 1)
		C.nlopt_set_ftol_rel(n.lopt, (C.double)(n.locFtolRel))
		C.nlopt_set_ftol_abs(n.lopt, (C.double)(n.locFtolAbs))
		C.nlopt_set_xtol_rel(n.lopt, (C.double)(n.locXtolRel))
		C.nlopt_set_local_optimizer(n.gopt, n.lopt)
	}

	nID := registerObject(n)
	defer unregisterObject(nID)
	C.nlopt_set_max_objective(n.gopt, (C.nlopt_func)(unsafe.Pointer(C.callback_adaptor)), unsafe.Pointer(&nID))

	lb := make([]C.double, len(n.parameters))
	ub := make([]C.double, len(n.parameters))
	x := make([]C.double, len(n.parameters))
	for i, par := range n.parameters {
		val := par.Get()
		nMin := par.GetMin()
		nMax := par.GetMax()
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
		n.ftolRel, n.ftolAbs, n.xtolRel, iterations)
	C.nlopt_set_ftol_rel(n.gopt, (C.double)(n.ftolRel))
	C.nlopt_set_ftol_abs(n.gopt, (C.double)(n.ftolAbs))
	C.nlopt_set_xtol_rel(n.gopt, (C.double)(n.xtolRel))
	C.nlopt_set_maxeval(n.gopt, (C.int)(iterations))

	var maxf C.double

	if n.seed > 0 {
		C.nlopt_srand((C.ulong)(n.seed))
	}

	n.PrintHeader()
	n.optRes = C.nlopt_optimize(n.gopt, (*C.double)(unsafe.Pointer(&x[0])), &maxf)
	if n.optRes < 0 {
		log.Fatalf("nlopt failed with code: %v (%v)", n.optRes, returnStatus[n.optRes])
	} else {
		log.Infof("nlopt success with code: %d (%v)", n.optRes, returnStatus[n.optRes])
	}

	for i, par := range n.parameters {
		par.Set((float64)(x[i]))
	}

	n.maxL = (float64)(maxf)
	n.maxLPar = n.parameters.Values(n.maxLPar)
}

// Summary returns optimization summary (i.e. success/error, etc).
func (n *NLOPT) Summary() interface{} {
	s := n.BaseOptimizer.Summary().(baseOptimizerSummary)
	s.Status = struct {
		Code       C.nlopt_result `json:"code"`
		CodeString string         `json:"codeString"`
	}{
		n.optRes,
		returnStatus[n.optRes],
	}
	return s
}

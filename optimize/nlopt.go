package optimize

// #cgo LDFLAGS: -lnlopt
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

type NLOPT struct {
	BaseOptimizer
	parameters FloatParameters
	Optimizable
	gopt C.nlopt_opt
	//lopt C.nlopt_opt
	dH          float64
	stop        bool
	finishLBFGS bool
}

func NewNLOPT() (nlopt *NLOPT) {
	nlopt = &NLOPT{
		BaseOptimizer: BaseOptimizer{
			repPeriod: 10,
		},
		dH:          1e-6,
		finishLBFGS: true,
	}
	return
}

func (n *NLOPT) SetOptimizable(opt Optimizable) {
	n.Optimizable = opt
	n.parameters = opt.GetFloatParameters()
}

func (n *NLOPT) Run(iterations int) {
	n.PrintHeader(n.parameters)

	n.gopt = C.nlopt_create(C.NLOPT_GN_CRS2_LM, (C.uint)(len(n.parameters)))
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

	C.nlopt_set_xtol_rel(n.gopt, 1e-2)
	C.nlopt_set_ftol_rel(n.gopt, 1e-4)

	var maxf C.double

	res := C.nlopt_optimize(n.gopt, (*C.double)(unsafe.Pointer(&x[0])), &maxf)
	if res < 0 {
		log.Fatalf("nlopt failed with code: %v", res)
	}

	for i, par := range n.parameters {
		par.Set((float64)(x[i]))
	}

	log.Info("Using LBFGSB to imrpove result")
	opt := NewLBFGSB()
	opt.SetOptimizable(n.Optimizable)
	opt.sig = n.sig
	opt.i = n.i
	opt.suppressHeader = true
	opt.Run(iterations)

	// if !n.Quiet {
	// 	log.Info("Finished NLOPT")
	// 	log.Noticef("Maximum likelihood: %v", maxf)
	// 	log.Infof("Parameter  names: %v", n.parameters.NamesString())
	// 	log.Infof("Parameter values: %v", n.GetMaxLParameters())
	// }
	// n.PrintFinal(n.parameters)
}

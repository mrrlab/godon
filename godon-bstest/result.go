package main

import (
	"errors"
)

type optimizer struct {
	// MaxLnL is the maximum log likelihood.
	MaxLnL float64
	// MaxLParameters is the maximum likelihood parameter values.
	MaxLParameters map[string]float64
}

// result stores the results of single godon execution.
type result struct {
	Optimizers []optimizer
	Model      interface{}
}

func (r *result) lastOptimizer() (*optimizer, error) {
	if len(r.Optimizers) == 0 {
		return nil, errors.New("optimizers array is empty")
	}
	return &r.Optimizers[len(r.Optimizers)-1], nil
}

// MustGetLikelihood returns likelihood or exits with an error.
func (r *result) MustGetLikelihood() float64 {
	res, err := r.GetLikelihood()
	if err != nil {
		log.Fatal("Cannot get likelihood:", err)
	}
	return res
}

// GetLikelihood gets maximum likelihood from the result.
func (r *result) GetLikelihood() (float64, error) {
	opt, err := r.lastOptimizer()
	if err != nil {
		return 0, err
	}
	return opt.MaxLnL, nil
}

// MustGetMaxLParameters returns maximum likelihood parameters or exits with an error.
func (r *result) MustGetMaxLParameters() map[string]float64 {
	res, err := r.GetMaxLParameters()
	if err != nil {
		log.Fatal("Cannot get maximum likelihood parameters:", err)
	}
	return res
}

// GetMaxLParameters gets maximum likelihood parameters from the result.
func (r *result) GetMaxLParameters() (map[string]float64, error) {
	opt, err := r.lastOptimizer()
	if err != nil {
		return nil, err
	}
	return opt.MaxLParameters, nil
}

// MustToHyp converts result to pointer to hyp, or exits in case of error.
func (r *result) MustToHyp() *hyp {
	res, err := r.ToHyp()
	if err != nil {
		log.Fatal("Cannot convert to hyp:", err)
	}
	return res
}

// ToHyp converts result to pointer to hyp.
func (r *result) ToHyp() (*hyp, error) {
	opt, err := r.lastOptimizer()
	if err != nil {
		return nil, err
	}
	return &hyp{
		MaxLnL:         opt.MaxLnL,
		MaxLParameters: opt.MaxLParameters,
		Final:          r.Model,
	}, nil

}

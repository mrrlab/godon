package main

type optimizer struct {
	// MaxLnL is the maximum log likelihood.
	MaxLnL float64
	// MaxLParameters is the maximum likelihood parameter values.
	MaxLParameters map[string]float64
}

// result stores the results of single godon execution.
type result struct {
	Optimizer optimizer
	Model     interface{}
}

// GetLikelihood gets maximum likelihood from the result.
func (r *result) GetLikelihood() float64 {
	return r.Optimizer.MaxLnL
}

// GetMaxLParameters gets maximum likelihood parameters from the result.
func (r *result) GetMaxLParameters() map[string]float64 {
	return r.Optimizer.MaxLParameters
}

// ToHyp converts result to pointer to hyp.
func (r *result) ToHyp() *hyp {
	return &hyp{
		MaxLnL:         r.Optimizer.MaxLnL,
		MaxLParameters: r.Optimizer.MaxLParameters,
		Final:          r.Model,
	}
}

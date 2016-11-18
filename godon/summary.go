package main

import "bitbucket.org/Davydov/godon/optimize"

type CallSummary struct {
	// Version stores godon version.
	Version string `json:"version"`
	// CommandLine is an array storing binary name and all command-line parameters.
	CommandLine []string `json:"commandLine"`
	// Seed is the seed used for random number generation initialization.
	Seed int64 `json:"seed"`
	// NThreads is the number of processes used.
	NThreads int `json:"nThreads"`
	// Time is the computations time in seconds.
	TotalTime float64 `json:"time"`
}

// OptimizationSummary is storing godon run summary information.
type OptimizationSummary struct {
	// StartingTree is the starting tree after unrooting.
	StartingTree string `json:"startingTree"`
	// FinalTree is the tree after branch length optimization (if performend).
	FinalTree string `json:"finalTree,omitempty"`
	// FullLnL is the full (non-aggregated) likelihood, it's only computed if specified (-printfull).
	FullLnL float64 `json:"fullLnL,omitempty"`
	// Time is the computations time in seconds.
	Time float64 `json:"optimizationTime"`
	// Model is the model summary, including BEB and NEB if available.
	Model interface{} `json:"model,omitempty"`
	// Optimizers is an array of summary of all optimizers used.
	Optimizer optimize.Summary `json:"optimizer"`
	// Hypothesis is H0 or H1
	Hypothesis string `json:"hypothesis,omitempty"`
}

// HypTestSummary is storing summary information for hypothesis test.
type HypTestSummary struct {
	// H0 is the result of H0 run.
	H0 HypSummary
	// H1 is the result of H1 run.
	H1 HypSummary
	// Runs stores all the runs
	Runs []OptimizationSummary `json:"runs"`
}

// Hyp summary stores information on one hypothesis
type HypSummary struct {
	// MaxLnL is the maximum log likelihood.
	MaxLnL float64 `json:"maxLnL"`
	// MaxLParameters is the maximum likelihood parameter values.
	MaxLParameters map[string]float64 `json:"maxLParameters"`
	// Model is the model summary, including BEB and NEB if available.
	Final interface{} `json:"final,omitempty"`
}

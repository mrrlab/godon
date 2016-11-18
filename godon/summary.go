package main

// OptimizationSummary is storing godon run summary information.
type OptimizationSummary struct {
	// Version stores godon version.
	Version string `json:"version"`
	// CommandLine is an array storing binary name and all command-line parameters.
	CommandLine []string `json:"commandLine"`
	// Seed is the seed used for random number generation initialization.
	Seed int64 `json:"seed"`
	// NThreads is the number of processes used.
	NThreads int `json:"nThreads"`
	// StartingTree is the starting tree after unrooting.
	StartingTree string `json:"startingTree"`
	// FinalTree is the tree after branch length optimization (if performend).
	FinalTree string `json:"finalTree,omitempty"`
	// FullLnL is the full (non-aggregated) likelihood, it's only computed if specified (-printfull).
	FullLnL float64 `json:"fullLnL,omitempty"`
	// Time is the computations time in seconds.
	Time float64 `json:"time"`
	// Model is the model summary, including BEB and NEB if available.
	Model interface{} `json:"model,omitempty"`
	// Optimizers is an array of summary of all optimizers used.
	Optimizer interface{} `json:"optimizer"`
}

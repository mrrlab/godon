package main

// Summary is storing godon summary information.
type Summary struct {
	// Version stores godon version.
	Version string `json:"version"`
	// CommandLine is an array storing binary name and all command-line parameters.
	CommandLine []string `json:"commandLine"`
	// Seed is the seed used for random number generation initialization.
	Seed int64 `json:"seed"`
	// NThreads is the number of processes used.
	NThreads int `json:"nThreads"`
	// StartTree is the starting tree after unrooting.
	StartTree string `json:"startTree"`
	// EndTree is the tree after branch length optimization (if performend).
	EndTree string `json:"endTree,omitempty"`
	// FullLnL is the full (non-aggregated) likelihood, it's only computed if specified (-printfull).
	FullLnL float64 `json:"fullLnL,omitempty"`
	// Time is the computations time in seconds.
	Time float64 `json:"time"`
	// Model is the model summary, including BEB and NEB if available.
	Model interface{} `json:"model,omitempty"`
	// Optimizers is an array of summary of all optimizers used.
	Optimizers []interface{} `json:"optimizers"`
}

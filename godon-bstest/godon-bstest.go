/*

Godon-bstest is a wrapper for godon. It runs godon multiple times to
perform a test, such as branch-site test for positive selection.

*/
package main

import (
	"encoding/json"
	"flag"
	"io/ioutil"
	"os"
	"path"
	"time"

	"github.com/op/go-logging"
)

// defaultSignificanceThreshold is the minimal D-statistics value, for which
// beb should be computed (qchisq(0.9, df=1)).
const defaultSignificanceThreshold = 2.705543

// defaultThreshold is the minimal likelihood improvement to
// consider an improvement for subsequent rerunning
const defaultUpdateThreshold = 1e-3

// setting up logging
var formatter = logging.MustStringFormatter(`%{message}`)
var log = logging.MustGetLogger("godon-bstest")

// dir is the temporary directory
var dir string

func saveJSON() {
	if *jsonF != "" {
		j, err := json.Marshal(sum)
		if err != nil {
			log.Error(err)
		} else {
			log.Debug(string(j))
			f, err := os.Create(*jsonF)
			if err != nil {
				log.Error("Error creating json output file:", err)
			} else {
				f.Write(j)
			}
			f.Close()
		}
	}
}

type hyp struct {
	// MaxLnL is the maximum log likelihood.
	MaxLnL float64 `json:"maxLnL"`
	// MaxLParameters is the maximum likelihood parameter values.
	MaxLParameters map[string]float64 `json:"maxLParameters"`
	// Model is the model summary, including BEB and NEB if available.
	Final interface{} `json:"final,omitempty"`
}

// summary is storing summary information.
type summary struct {
	// H0 is the result of H0 run.
	H0 *hyp
	// H0 is the result of H0 run.
	H1 *hyp
	// Runs stores all the runs
	Runs []interface{} `json:"runs"`
	// Time is the computations time in seconds.
	Time float64 `json:"time"`
}

// sum contains the run summary
var sum summary

// program parameters
var binary = flag.String("binary", "godon", "binary name or full path to godon")
var debug = flag.Bool("debug", false, "enable debug mode")
var jsonF = flag.String("json", "", "write json output to a file")
var sthr = flag.Float64("sthr", defaultSignificanceThreshold,
	"LRT siginficance threshold, for H0 rerun and posterior computations")
var uthr = flag.Float64("uthr", defaultUpdateThreshold, "likelihood improvement threshold for rerunning H0")
var quick = flag.Bool("quick", false, "only prevent negative LRT statistics")

func main() {
	startTime := time.Now()

	logging.SetFormatter(formatter)
	backend := logging.NewLogBackend(os.Stderr, "", 0)
	logging.SetBackend(backend)

	flag.Parse()

	if !*debug {
		logging.SetLevel(logging.INFO, "godon-bstest")
	}

	if len(flag.Args()) < 2 {
		log.Fatal("you should specify tree and alignment")
		return
	}

	var err error

	// create a temporary dir for files
	dir, err = ioutil.TempDir("", path.Base(flag.Arg(0)))
	if err != nil {
		log.Fatal(err)
	}
	defer os.RemoveAll(dir)

	res0 := mustRun(false, nil, false, true)
	res1 := mustRun(true, nil, false, true)

	var l0, l1 float64
	l0 = res0.GetLikelihood()
	l1 = res1.GetLikelihood()

	for updated := true; updated; {
		// stop if nothing has been updated.
		// this loop normally should run for a single time.
		updated = false
		// if l1 < l0, rerun H1 starting from H0.
		if l1 < l0 {
			updated = true
			h0par := res0.GetMaxLParameters()
			res1 = mustRun(true, h0par, false, !*quick)
			l1 = res1.GetLikelihood()
		}

		// if significant (D>thr), rerun H0 starting from H1
		if 2*(l1-l0) > *sthr && !*quick {
			h1par := res1.GetMaxLParameters()
			res0Alt := mustRun(false, h1par, false, true)
			l0Alt := res0Alt.GetLikelihood()
			if l0Alt > l0 {
				if l0Alt-l0 > *uthr {
					updated = true
				}
				res0 = res0Alt
				l0 = l0Alt
			}
		}
	}

	// final BEB & NEB computation
	if 2*(l1-l0) > *sthr {
		h1par := res1.GetMaxLParameters()
		res1 = mustRun(true, h1par, true, false)
	}

	sum.H0 = res0.ToHyp()
	sum.H1 = res1.ToHyp()

	log.Infof("lnL0=%f, lnL1=%f", sum.H0.MaxLnL, sum.H1.MaxLnL)
	endTime := time.Now()

	deltaT := endTime.Sub(startTime)
	log.Noticef("Running time: %v", deltaT)

	sum.Time = deltaT.Seconds()

	saveJSON()
}

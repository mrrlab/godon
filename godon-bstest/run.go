package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"path"
)

// counter is used to make unique filenames
var counter = 1

// constructArgs constructs arguments for running godon.
func constructArgs(h1 bool, jsonF string, start string, final bool) (args []string) {
	log.Debugf("output json to: %v", jsonF)
	args = make([]string, 0, len(flag.Args())+10)

	args = append(args, flag.Args()[2:]...)
	args = append(args, "-json", jsonF)
	if !h1 {
		args = append(args, "-fixw")
	}
	if len(start) > 0 {
		args = append(args, "-start", start)
	}
	// this should be after args, we are overriding method
	if !final {
		args = append(args, "-nofinal")
	} else {
		args = append(args, "-method", "none")
	}
	args = append(args, flag.Args()[:2]...)
	log.Debugf("%v (%v)", args, len(args))
	return
}

// mustRun runs optimization and exits with if an error is encountered.
func mustRun(h1 bool, start map[string]float64, final bool) *result {
	res, err := run(h1, start, final)
	if err != nil {
		log.Fatal("Error running godon:", err)
	}
	return res
}

// run runs optimization for specific parameters.
func run(h1 bool, start map[string]float64, final bool) (*result, error) {
	hyp := 0
	if h1 {
		hyp = 1
	}

	jBaseName := fmt.Sprintf("%03d_H%d", counter, hyp)
	counter++
	jFullName := path.Join(dir, jBaseName+".json")

	startFileName := ""

	if len(start) > 0 {
		start := copyMap(start)
		if h1 {
			if _, ok := start["omega2"]; !ok {
				start["omega2"] = 1
			}
		} else {
			delete(start, "omega2")
		}
		startFileName = path.Join(dir, jBaseName+"_start.json")
		j, err := json.Marshal(start)
		if err != nil {
			return nil, err
		}
		f, err := os.Create(startFileName)
		if err != nil {
			return nil, err
		}

		f.Write(j)
		f.Close()
	}

	args := constructArgs(h1, jFullName, startFileName, final)
	cmd := exec.Command(*binary, args...)
	cmd.Stderr = os.Stderr
	err := cmd.Run()
	if err != nil {
		return nil, err
	}
	fileContents, err := ioutil.ReadFile(jFullName)
	if err != nil {
		return nil, err
	}

	// save full results
	var fullRes interface{}
	err = json.Unmarshal(fileContents, &fullRes)
	if err != nil {
		return nil, err
	}
	sum.Runs = append(sum.Runs, fullRes)

	var res result
	err = json.Unmarshal(fileContents, &res)
	if err != nil {
		return nil, err
	}

	return &res, nil
}

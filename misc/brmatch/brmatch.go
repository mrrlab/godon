// Brmatch prints branch label (*label) correspondance between two
// topologically identical trees.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"bitbucket.org/Davydov/godon/tree"
)

func main() {
	flag.Parse()

	if flag.NArg() < 2 {
		log.Fatal("Please provide two files")
	}
	in1file, err := os.Open(flag.Arg(0))
	if err != nil {
		log.Fatal(err)
	}
	defer in1file.Close()

	in2file, err := os.Open(flag.Arg(1))
	if err != nil {
		log.Fatal(err)
	}
	defer in2file.Close()

	t1, err := tree.ParseNewick(in1file)
	if err != nil {
		log.Fatal(err)
	}

	t2, err := tree.ParseNewick(in2file)
	if err != nil {
		log.Fatal(err)
	}

	nodes1 := t1.NodeIDArray()
	nodes2 := t2.NodeIDArray()
	if len(nodes1) != len(nodes2) {
		log.Fatal("Two trees have different number of nodes")
	}

	for i := range nodes1 {
		nm1 := nodes1[i].Name
		if nm1 != "" && nm1[0] == '*' {
			nm2 := nodes2[i].Name
			if nm2 != "" && nm2[0] == '*' {
				fmt.Println(nm1[1:], nm2[1:])
			}
		}
	}
}

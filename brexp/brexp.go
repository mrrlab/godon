package main

import (
	"fmt"
	"log"
	"os"

	"bitbucket.com/Davydov/golh/tree"
)

func main() {
	t, err := tree.ParseNewick(os.Stdin)
	if err != nil {
		log.Fatal(err)
	}
	for _, node := range t.Nodes() {
		fmt.Printf("br%d=%f\n", node.Id, node.BranchLength)
	}
}

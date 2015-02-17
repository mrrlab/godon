package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"regexp"

	"bitbucket.com/Davydov/golh/tree"
)

func SelectomeBrString(t *tree.Tree) (s string) {
	s = SelectomeNodeStringBr(t.Node)
	n := 1
	repfunc := func(sym string) string {
		s := fmt.Sprintf("*%d", n)
		n++
		return s
	}
	es := regexp.MustCompile(`\*|!`)
	s = es.ReplaceAllStringFunc(s, repfunc)
	return

}

func SelectomeNodeStringBr(node *tree.Node) (s string) {
	if node.IsTerminal() {
		return fmt.Sprintf("%s", node.Name)
	}
	s += "("
	for i, child := range node.ChildNodes() {
		s += SelectomeNodeStringBr(child)
		if i != len(node.ChildNodes())-1 {
			s += ","
		}
	}
	s += ")"
	if node.IsRoot() {
		s += ";"
	} else {
		s += "*"
	}
	return s
}

func main() {
	infilename := flag.String("in", "", "input filename")
	mode := flag.String("mode", "brlen", "program mode (default: brlen)")

	flag.Parse()

	var infile *os.File
	var err error

	if *infilename == "" {
		infile = os.Stdin
	} else {
		infile, err = os.Open(*infilename)
		if err != nil {
			log.Fatal(err)
		}
		defer infile.Close()
	}

	t, err := tree.ParseNewick(infile)
	if err != nil {
		log.Fatal(err)
	}
	switch *mode {
	case "brlen":
		for _, node := range t.NodeIdArray() {
			if node != nil {
				fmt.Printf("br%d=%f\n", node.Id, node.BranchLength)
			}
		}
	case "brtree":
		fmt.Println(t.BrString())
	case "selectome":
		fmt.Println(SelectomeBrString(t))
	default:
		log.Fatal("Unknown mode")
	}
}

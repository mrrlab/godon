package main

import (
	"bio"
	"fmt"
	"math"
	"os"
	"runtime"
	"tree"
)

var (
	alphabet = [...]byte{'A', 'T', 'G', 'C'}
	nm2id    = make(map[string]int)
)

func L(ali bio.Sequences, t *tree.Tree) (lnL float64) {
	ch := make(chan float64, len(ali[0].Sequence))
	for i, _ := range ali[0].Sequence {
		go func(i int) {
			//fmt.Println(i);
			ch <- subL(ali, t, i)
		}(i)
	}

	for _, _ = range ali[0].Sequence {
		lnL += <-ch
	}
	close(ch)
	return
}

func subL(ali bio.Sequences, t *tree.Tree, i int) float64 {
	//fmt.Println(i)
	res := 0.0
	plh := make(map[*tree.Tree]map[byte]float64)
	for node := range t.Nodes() {
		plh[node] = make(map[byte]float64)
	}

	nodes := make(chan *tree.Tree, len(ali))
	for node := range t.Terminals() {
		for _, l := range alphabet {
			if l == ali[nm2id[node.Name]].Sequence[i] {
				plh[node][l] = 1
			} else {
				plh[node][l] = 0
			}
		}
		nodes <- node.Parent
	}

NodeLoop:
	for node := range nodes {
		for child := range node.ChildNodes() {
			if len(plh[child]) == 0 {
				nodes <- node
				continue NodeLoop
			}
		}
		if len(plh[node]) > 0 {
			continue NodeLoop
		}
		for _, l1 := range alphabet {
			l := 1.0
			for child := range node.ChildNodes() {
				s := 0.0
				for _, l2 := range alphabet {
					var x float64
					if l1 == l2 {
						x = 1./4 + 3./4*math.Exp(-4.*child.BranchLength/3)
					} else {
						x = 1./4 - 1./4*math.Exp(-4.*child.BranchLength/3)
					}
					s += x * plh[child][l2]
				}
				l *= s
			}
			plh[node][l1] = l
		}
		nodes <- node.Parent
		if node.IsRoot() {
			for _, l := range alphabet {
				res += 0.25 * plh[node][l]
			}
			break NodeLoop
		}

	}
	close(nodes)
	return math.Log(res)
}

func main() {
	if len(os.Args) < 3 {
		fmt.Println("specify files")
		return
	}
	t, err := tree.ParseNewick(os.Args[1])
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Println(t.FullString())
	ali, err := bio.ParseFasta(os.Args[2])
	if err != nil {
		fmt.Println(err)
		return
	}

	runtime.GOMAXPROCS(2)
	fmt.Println("Will use ", runtime.GOMAXPROCS(0), "CPUs")

	for i, s := range ali {
		nm2id[s.Name] = i
	}

	fmt.Println(nm2id)

	for i := 0; i < 100; i++ {
		fmt.Println("lnL=", L(ali, t))
	}

	//var lnL float64

	//fmt.Println(lnL)
}

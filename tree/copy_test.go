package tree

import (
	"bytes"
	"testing"
)

const (
	tree1 = "((((a001:0.242690,a002:0.268555)#1:0.073424,a003:0.252510):0.198740,((((((a004:0.001000,a005:0.014869):0.045007,a006:0.050606):0.056908,a007:0.166439):0.023217,a008:0.094788):0.429852,a009:0.558116):0.130317,(a010:0.009332,a011:0.024271):0.315124):0.217376):0.464470,a012:0.144369):0.0;"
)

func TestCopy1(tst *testing.T) {
	t, err := ParseNewick(bytes.NewBufferString(tree1))
	if err != nil {
		tst.Error("Error parsing tree", err)
	}
	t1 := t.Copy()
	t2 := t1.Copy()

	t.ClearCache()
	t1.ClearCache()

	tNodes := t.Nodes()
	t1Nodes := t1.Nodes()
	t2Nodes := t2.Nodes()

	if len(tNodes) != len(t1Nodes) {
		tst.Error("node length differ between t and t1")
	}
	if len(t1Nodes) != len(t2Nodes) {
		tst.Error("node length differ between t1 and t2")
	}

	for i := 0; i < len(tNodes); i++ {
		if tNodes[i] == t1Nodes[i] ||
			t1Nodes[i] == t2Nodes[i] {
			tst.Error("node pointers match between trees")
		}
		if tNodes[i].BranchLength != t1Nodes[i].BranchLength ||
			t1Nodes[i].BranchLength != t2Nodes[i].BranchLength {
			tst.Error("node length differ")
		}
		if tNodes[i].Name != t1Nodes[i].Name ||
			t1Nodes[i].Name != t2Nodes[i].Name {
			tst.Error("node name differ")
		}
		if tNodes[i].Class != t1Nodes[i].Class ||
			t1Nodes[i].Class != t2Nodes[i].Class {
			tst.Error("node name differ")
		}
	}

	for _, node := range t1.Nodes() {
		node.BranchLength = 2
	}

	for i := 0; i < len(tNodes); i++ {
		if t.Nodes()[i].BranchLength == t1.Nodes()[i].BranchLength {
			tst.Error("node length still match after change")
		}
	}

	for _, node := range t2.Nodes() {
		node.BranchLength = 0.5
	}

	for i := 0; i < len(tNodes); i++ {
		if t1.Nodes()[i].BranchLength <= t2.Nodes()[i].BranchLength {
			tst.Error("node length is wrong")
		}
	}
}

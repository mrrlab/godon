package tree

import (
	"bytes"
	"testing"
)

const (
	tree2 = "((a:1,b:2)#1:3,c:1):0;"
	tree3 = "(c#1:1,(a:1,b:2):3):0;"
)

func TestUnroot1(tst *testing.T) {
	t, err := ParseNewick(bytes.NewBufferString(tree2))
	if err != nil {
		tst.Error("Error parsing tree", err)
	}
	tst.Log("Got tree:", t.ClassString())

	// Unroot
	b, err := t.Unroot()
	if err != nil {
		tst.Error("Error unrooting tree", err)
	}
	tst.Log("Unrooted:", t.ClassString())
	tst.Log("b =", b)
	if t.String() != "(a:1.000000,b:2.000000,c:4.000000):0.000000;" {
		tst.Error("Error unrooting tree, got:", t)
	}

	err = t.Root(b)
	if err != nil {
		tst.Error("Error rooting tree", err)
	}

	tst.Log("Rooted:", t.ClassString())
	if t.String() != "((a:1.000000,b:2.000000):2.000000,c:2.000000):0.000000;" {
		tst.Error("Error unrooting tree, got:", t)
	}
}

func TestUnroot2(tst *testing.T) {
	t, err := ParseNewick(bytes.NewBufferString(tree3))
	if err != nil {
		tst.Error("Error parsing tree", err)
	}
	tst.Log("Got tree:", t.ClassString())

	// Unroot
	b, err := t.Unroot()
	if err != nil {
		tst.Error("Error unrooting tree", err)
	}
	tst.Log("Unrooted:", t.ClassString())
	tst.Log("b =", b)
	if t.String() != "(c:4.000000,a:1.000000,b:2.000000):0.000000;" {
		tst.Error("Error unrooting tree, got:", t)
	}

	err = t.Root(b)
	if err != nil {
		tst.Error("Error rooting tree", err)
	}

	tst.Log("Rooted:", t.ClassString())
	if t.String() != "(c:2.000000,(a:1.000000,b:2.000000):2.000000):0.000000;" {
		tst.Error("Error unrooting tree, got:", t)
	}
}

func TestDouble1(tst *testing.T) {
	t, err := ParseNewick(bytes.NewBufferString(tree1))
	if err != nil {
		tst.Error("Error parsing tree", err)
	}
	tst.Log("Got tree:", t)

	// Unroot
	b, err := t.Unroot()
	if err != nil {
		tst.Error("Error unrooting tree", err)
	}
	tunr := t.String()
	tst.Log("Unrooted:", t)

	err = t.Root(b)
	if err != nil {
		tst.Error("Error rooting tree", err)
	}
	tr := t.String()
	tst.Log("Rooted:", t)

	b, err = t.Unroot()
	if err != nil {
		tst.Error("Error unrooting tree", err)
	}
	tst.Log("Unrooted:", t)
	if t.String() != tunr {
		tst.Error("Multiple rooting/unrooting fail")
	}

	err = t.Root(b)
	if err != nil {
		tst.Error("Error rooting tree", err)
	}
	tst.Log("Rooted:", t)
	if t.String() != tr {
		tst.Error("Multiple rooting/unrooting fail")
	}
}

//TODO tests for class propagation and two-leaves error

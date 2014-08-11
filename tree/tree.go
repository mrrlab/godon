package tree

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"
	"unicode"
	"unicode/utf8"
)

var nodeId int

type Mode int

const (
	NORMAL Mode = iota
	LENGTH
	CLASS
)

type Tree struct {
	Name         string
	BranchLength float64
	Parent       *Tree
	childNodes   []*Tree
	Id           int
	Class        int
}

func NewTree(parent *Tree) (tree *Tree) {
	tree = &Tree{Parent: parent, Id: nodeId}
	nodeId++
	return
}

func (tree *Tree) AddChild(subTree *Tree) {
	/*if tree.childNodes == nil{
		tree.childNodes = make([]*Tree, 0, 2)
	}*/
	subTree.Parent = tree
	tree.childNodes = append(tree.childNodes, subTree)
}

func (tree *Tree) String() (s string) {
	s = "<"
	if tree.Parent == nil {
		s += "root, "
	}
	if tree.Name != "" {
		s += "name=" + tree.Name + ", "
	}
	s += fmt.Sprintf("Id=%v, BranchLength=%v", tree.Id, tree.BranchLength)
	if tree.Class != 0 {
		s += fmt.Sprintf(", Class=%v", tree.Class)
	}
	s += ">"
	return
}

func (tree *Tree) FullString() string {
	return strings.TrimSpace(tree.prefixString(""))
}

func (tree *Tree) prefixString(prefix string) (s string) {
	s = prefix + tree.String() + "\n"
	for _, node := range tree.childNodes {
		s += node.prefixString(prefix + "    ")
	}
	return
}

func (tree *Tree) ChildNodes() <-chan *Tree {
	ch := make(chan *Tree)
	go func() {
		for _, node := range tree.childNodes {
			ch <- node
		}
		close(ch)
	}()
	return ch
}

func (tree *Tree) Walk(ch chan *Tree, filter func(*Tree) bool) {
	if filter == nil || filter(tree) {
		ch <- tree
	}
	for _, node := range tree.childNodes {
		node.Walk(ch, filter)
	}
}

func (tree *Tree) Walker(filter func(*Tree) bool) <-chan *Tree {
	ch := make(chan *Tree)
	go func() {
		tree.Walk(ch, filter)
		close(ch)
	}()
	return ch
}

func (tree *Tree) Nodes() <-chan *Tree {
	return tree.Walker(nil)
}

func (tree *Tree) NNodes() (size int) {
	for _, node := range tree.childNodes {
		size += node.NNodes()
	}
	return size + 1
}
func (tree *Tree) Terminals() <-chan *Tree {
	return tree.Walker(func(t *Tree) bool {
		if len(t.childNodes) == 0 {
			return true
		}
		return false
	})
}

func (tree *Tree) NonTerminals() <-chan *Tree {
	return tree.Walker(func(t *Tree) bool {
		if len(t.childNodes) > 0 {
			return true
		}
		return false
	})
}

func (tree *Tree) ClassNodes(class int) <-chan *Tree {
	return tree.Walker(func(t *Tree) bool {
		if t.Class == class {
			return true
		}
		return false
	})
}

func (tree *Tree) IsRoot() bool {
	return tree.Parent == nil
}

func IsSpecial(c rune) bool {
	switch c {
	case '(', ')', ':', '#', ';', ',':
		return true
	}
	return false

}
func NewickSplit(data []byte, atEOF bool) (advance int, token []byte, err error) {
	start := 0
	// Skip leading spaces; and return 1-char tokens.
	for width := 0; start < len(data); start += width {
		var r rune
		r, width = utf8.DecodeRune(data[start:])
		if IsSpecial(r) {
			return start + width, data[start : start+width], nil
		}
		if !unicode.IsSpace(r) {
			break
		}
	}
	if atEOF && len(data) == 0 {
		return 0, nil, nil
	}

	// Scan until space or special character.
	for width, i := 0, start; i < len(data); i += width {
		var r rune
		r, width = utf8.DecodeRune(data[i:])
		if unicode.IsSpace(r) || IsSpecial(r) {
			return i, data[start:i], nil
		}
	}
	// If we're at EOF, we have a final, non-empty, non-terminated word. Return it.
	if atEOF && len(data) > start {
		return len(data), data[start:], nil
	}
	// Request more data.
	return 0, nil, nil
}

func ParseNewick(rd io.Reader) (tree *Tree, err error) {
	scanner := bufio.NewScanner(rd)

	scanner.Split(NewickSplit)

	nodeId = 0

	tree = NewTree(nil)
	mode := NORMAL

	for scanner.Scan() {
		text := scanner.Text()
		switch text {
		case "(":
			subTree := NewTree(nil)
			if tree != nil {
				tree.AddChild(subTree)
			}
			tree = subTree

		case ",":
			if tree.Parent == nil {
				return nil, errors.New("top level comma mismatch")
			}
			subTree := NewTree(nil)
			tree.Parent.AddChild(subTree)
			tree = subTree

		case ")":
			if tree.Parent == nil {
				return nil, errors.New("brackets mismatch")
			}
			tree = tree.Parent
		case "#":
			mode = CLASS
		case ":":
			mode = LENGTH
		case ";":
			return
		default:
			switch mode {
			case LENGTH:
				l, err := strconv.ParseFloat(text, 64)
				if err != nil {
					return nil, err
				}
				tree.BranchLength = l
				mode = NORMAL
			case CLASS:
				cl, err := strconv.ParseInt(text, 0, 0)
				if err != nil {
					return nil, err
				}
				tree.Class = int(cl)
				mode = NORMAL
			default:
				tree.Name = text
			}
		}
	}

	return

}

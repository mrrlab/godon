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

type Mode int

const (
	NORMAL Mode = iota
	LENGTH
	CLASS
)

type Tree struct {
	*Node
	nNodes    int
	nodes     []*Node
	nodeOrder []*Node
}

func (tree *Tree) ClearCache() {
	tree.nNodes = 0
	// can we use something like append here? do we need to bother?
	tree.nodes = nil
	tree.nodeOrder = nil
}

func (tree *Tree) NNodes() int {
	if tree.nNodes == 0 {
		tree.nNodes = tree.NSubNodes()
	}
	return tree.nNodes
}

func (tree *Tree) Nodes() []*Node {
	if tree.nodes == nil {
		tree.nodes = make([]*Node, tree.NNodes())
		i := 0
		for node := range tree.Walker(nil) {
			tree.nodes[node.Id] = node
			i++
		}
	}
	return tree.nodes
}

func (tree *Tree) Terminals() <-chan *Node {
	return tree.Walker(func(n *Node) bool {
		if len(n.childNodes) == 0 {
			return true
		}
		return false
	})
}

func (tree *Tree) NonTerminals() <-chan *Node {
	return tree.Walker(func(node *Node) bool {
		return node.IsTerminal()
	})
}

func (tree *Tree) ClassNodes(class int) <-chan *Node {
	return tree.Walker(func(node *Node) bool {
		if node.Class == class {
			return true
		}
		return false
	})
}

func (tree *Tree) NLeaves() (i int) {
	for _ = range tree.Walker(func(node *Node) bool {
		return node.IsTerminal()
	}) {
		i++
	}
	return
}

func (tree *Tree) Walker(filter func(*Node) bool) <-chan *Node {
	ch := make(chan *Node, tree.NNodes())
	tree.Walk(ch, filter)
	close(ch)
	return ch
}

// Copy creates independent copy of the tree.
func (tree *Tree) Copy() (newTree *Tree) {
	nNodes := tree.NNodes()
	newTree = &Tree{
		nNodes:    nNodes,
		nodes:     make([]*Node, nNodes),
		nodeOrder: make([]*Node, len(tree.NodeOrder())),
	}

	// Create node list.
	for i, node := range tree.Nodes() {
		if i != node.Id {
			panic("node id mismatch")
		}
		newTree.nodes[i] = node.Copy()
	}

	// Rewire node/parent connections.
	for i, node := range tree.Nodes() {
		newNode := newTree.nodes[i]
		for _, child := range node.childNodes {
			newChild := newTree.nodes[child.Id]
			newNode.AddChild(newChild)
		}
	}

	// Set nodeOrder
	for i, node := range tree.NodeOrder() {
		newTree.nodeOrder[i] = newTree.nodes[node.Id]
	}

	// Set root node.
	newTree.Node = newTree.nodes[0]

	return
}

func (tree *Tree) NodeOrder() []*Node {
	if tree.nodeOrder == nil {
		tree.nodeOrder = make([]*Node, 0, tree.NNodes())
		computed := make(map[*Node]bool, tree.NNodes())
		awaiting := make(chan *Node, tree.NNodes()*2)
		for node := range tree.Terminals() {
			computed[node] = true
			awaiting <- node.Parent
		}

		for node := range awaiting {
			if node == nil {
				break
			}
			if computed[node] {
				continue
			}
			allComputed := true
			for _, childNode := range node.ChildNodes() {
				if !computed[childNode] {
					allComputed = false
					break
				}
			}
			if !allComputed {
				awaiting <- node
			} else {
				tree.nodeOrder = append(tree.nodeOrder, node)
				computed[node] = true
				awaiting <- node.Parent
			}
		}
	}
	return tree.nodeOrder
}

type Node struct {
	Name         string
	BranchLength float64
	Parent       *Node
	childNodes   []*Node
	Id           int
	LeafId       int
	Class        int
}

func NewNode(parent *Node, nodeId int) (node *Node) {
	node = &Node{Parent: parent, Id: nodeId}
	return
}

// Copy creates copy of node with empty parent and children.
func (node *Node) Copy() *Node {
	return &Node{
		Name:         node.Name,
		BranchLength: node.BranchLength,
		childNodes:   make([]*Node, 0, len(node.childNodes)),
		Id:           node.Id,
		LeafId:       node.LeafId,
		Class:        node.Class,
	}
}

func (node *Node) AddChild(subNode *Node) {
	subNode.Parent = node
	node.childNodes = append(node.childNodes, subNode)
}

func (node *Node) StringBr() (s string) {
	if node.IsTerminal() {
		return fmt.Sprintf("%s#br%d", node.Name, node.Id)
	}
	s += "("
	for i, child := range node.childNodes {
		s += child.StringBr()
		if i != len(node.childNodes)-1 {
			s += ","
		}
	}
	s += fmt.Sprintf(")#br%d", node.Id)
	if node.IsRoot() {
		s += ";"
	}
	return s
}

func (node *Node) String() (s string) {
	if node.IsTerminal() {
		return fmt.Sprintf("%s:%0.6f", node.Name, node.BranchLength)
	}
	s += "("
	for i, child := range node.childNodes {
		s += child.String()
		if i != len(node.childNodes)-1 {
			s += ","
		}
	}
	s += fmt.Sprintf("):%0.6f", node.BranchLength)
	if node.IsRoot() {
		s += ";"
	}
	return s
}

func (node *Node) LongString() (s string) {
	s = "<"
	if node.Parent == nil {
		s += "root, "
	}
	if node.Name != "" {
		s += "name=" + node.Name + ", "
	}
	s += fmt.Sprintf("Id=%v, BranchLength=%v", node.Id, node.BranchLength)
	if node.IsTerminal() {
		s += fmt.Sprintf(", TipId=%v", node.LeafId)
	}
	if node.Class != 0 {
		s += fmt.Sprintf(", Class=%v", node.Class)
	}
	s += ">"
	return
}

func (node *Node) FullString() string {
	return strings.TrimSpace(node.prefixString(""))
}

func (node *Node) prefixString(prefix string) (s string) {
	s = prefix + node.LongString() + "\n"
	for _, node := range node.childNodes {
		s += node.prefixString(prefix + "    ")
	}
	return
}

func (node *Node) ChildNodes() []*Node {
	return node.childNodes
}

func (node *Node) Walk(ch chan *Node, filter func(*Node) bool) {
	if filter == nil || filter(node) {
		ch <- node
	}
	for _, node := range node.childNodes {
		node.Walk(ch, filter)
	}
}

func (node *Node) NSubNodes() (size int) {
	for _, node := range node.childNodes {
		size += node.NSubNodes()
	}
	return size + 1
}

func (node *Node) IsRoot() bool {
	return node.Parent == nil
}

func (node *Node) IsTerminal() bool {
	return len(node.childNodes) == 0
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

	nodeId := 0
	leafId := 0

	node := NewNode(nil, nodeId)
	tree = &Tree{Node: node}
	nodeId++

	mode := NORMAL

	for scanner.Scan() {
		text := scanner.Text()
		switch text {
		case "(":
			subNode := NewNode(nil, nodeId)
			nodeId++
			if node != nil {
				node.AddChild(subNode)
			}
			node = subNode

		case ",":
			if node.Parent == nil {
				return nil, errors.New("top level comma mismatch")
			}
			subNode := NewNode(nil, nodeId)
			nodeId++

			node.Parent.AddChild(subNode)
			node = subNode

		case ")":
			if node.Parent == nil {
				return nil, errors.New("brackets mismatch")
			}
			node = node.Parent
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
				node.BranchLength = l
				mode = NORMAL
			case CLASS:
				cl, err := strconv.ParseInt(text, 0, 0)
				if err != nil {
					return nil, err
				}
				node.Class = int(cl)
				mode = NORMAL
			default:
				node.LeafId = leafId
				leafId++
				node.Name = text
			}
		}
	}

	return

}

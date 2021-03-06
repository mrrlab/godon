// Package tree implements tree structure and newick parsing.
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

// Mode is a newick parsing mode variable.
type Mode int

// Newick parsing modes.
const (
	// Default parsing mode.
	NORMAL Mode = iota
	// Parsing branch length.
	LENGTH
	// Parsing node class.
	CLASS
)

// Tree is a tree structure.
type Tree struct {
	// Node is the root node.
	*Node
	nNodes      int
	maxNodeID   int
	nodeIDArray []*Node
	nodeOrder   []*Node
}

// ClearCache clears all precomputed parameters.
func (tree *Tree) ClearCache() {
	tree.nNodes = 0
	// can we use something like append here? do we need to bother?
	tree.nodeIDArray = nil
	tree.nodeOrder = nil
}

// NNodes returns number of nodes.
func (tree *Tree) NNodes() int {
	if tree.nNodes == 0 {
		tree.nNodes = tree.NSubNodes()
	}
	return tree.nNodes
}

// MaxNodeID returns maximum id of a node.
func (tree *Tree) MaxNodeID() (maxID int) {
	if tree.maxNodeID > 0 {
		return tree.maxNodeID
	}
	for node := range tree.Walker(nil) {
		if maxID < node.ID {
			maxID = node.ID
		}
	}
	tree.maxNodeID = maxID
	return
}

// NodeIDArray returns an array (slice). This array allows addressing
// nodes by their ids.
func (tree *Tree) NodeIDArray() []*Node {
	if tree.nodeIDArray == nil {
		tree.nodeIDArray = make([]*Node, tree.MaxNodeID()+1)
		i := 0
		for node := range tree.Walker(nil) {
			tree.nodeIDArray[node.ID] = node
			i++
		}
	}
	return tree.nodeIDArray
}

// Terminals returns a channel with all the terminal nodes of a tree.
func (tree *Tree) Terminals() <-chan *Node {
	return tree.Walker(func(n *Node) bool {
		return len(n.childNodes) == 0
	})
}

// NonTerminals returns a channel with all the non-terminal nodes of a
// tree.
func (tree *Tree) NonTerminals() <-chan *Node {
	return tree.Walker(func(node *Node) bool {
		return node.IsTerminal()
	})
}

// ClassNodes returns a channel with all the nodes of a given class.
func (tree *Tree) ClassNodes(class int) <-chan *Node {
	return tree.Walker(func(node *Node) bool {
		return node.Class == class
	})
}

// NLeaves returns number of leaves.
func (tree *Tree) NLeaves() (i int) {
	for range tree.Walker(func(node *Node) bool {
		return node.IsTerminal()
	}) {
		i++
	}
	return
}

// Walker returns a channel iterating over nodes. An optional function
// is provided to filter nodes.
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
		nNodes:      nNodes,
		nodeIDArray: make([]*Node, len(tree.NodeIDArray())),
		nodeOrder:   make([]*Node, len(tree.NodeOrder())),
	}

	// Create node list.
	for _, node := range tree.NodeIDArray() {
		if node != nil {
			newTree.nodeIDArray[node.ID] = node.Copy()
		}
	}

	// Rewire node/parent connections.
	for i, node := range tree.NodeIDArray() {
		if node == nil {
			continue
		}
		newNode := newTree.nodeIDArray[i]
		for _, child := range node.childNodes {
			newChild := newTree.nodeIDArray[child.ID]
			newNode.AddChild(newChild)
		}
	}

	// Set nodeOrder
	for i, node := range tree.NodeOrder() {
		newTree.nodeOrder[i] = newTree.nodeIDArray[node.ID]
	}

	// Set root node.
	newTree.Node = newTree.nodeIDArray[tree.Node.ID]

	return
}

// NodeOrder returns an array (slice) with non-terminal nodes in
// post-order.
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

// IsRooted returns true if the tree is rooted (i.e. bifurcation at
// the root).
func (tree *Tree) IsRooted() bool {
	return len(tree.Node.childNodes) == 2
}

// Unroot unroots the tree (i.e. replaces bifurcation at the root with
// a trifurcation).
func (tree *Tree) Unroot() (int, error) {
	if tree.NLeaves() < 3 {
		return 0, errors.New("Unrooting works only for trees with >= 2 leaves")
	}
	if len(tree.Node.childNodes) != 2 {
		// tree is probably unrooted already
		return 0, nil
	}

	reverse := false

	var newRoot, noRoot *Node
	switch {
	case !tree.Node.childNodes[0].IsTerminal():
		newRoot = tree.Node.childNodes[0]
		noRoot = tree.Node.childNodes[1]
	case !tree.Node.childNodes[1].IsTerminal():
		newRoot = tree.Node.childNodes[1]
		noRoot = tree.Node.childNodes[0]
		reverse = true
	default:
		// This probably shouldn't happen
		return 0, errors.New("Unexpected tree topology while unrooting")
	}

	if newRoot.Class != 0 && noRoot.Class != 0 && newRoot.Class != noRoot.Class {
		return 0, errors.New("Root branches have different non-zero classes")
	}

	tree.Node = newRoot
	noRoot.BranchLength += newRoot.BranchLength
	newRoot.BranchLength = 0
	if !reverse {
		newRoot.childNodes = append(newRoot.childNodes, noRoot)
	} else {
		newRoot.childNodes = append(newRoot.childNodes, nil)
		copy(newRoot.childNodes[1:], newRoot.childNodes)
		newRoot.childNodes[0] = noRoot
	}
	// rewire parents
	newRoot.Parent = nil
	noRoot.Parent = newRoot
	if newRoot.Class != 0 {
		noRoot.Class = newRoot.Class
		newRoot.Class = 0
	}

	// update cached parameters
	tree.nNodes--
	tree.nodeIDArray = nil
	tree.nodeOrder = nil

	return noRoot.ID, nil
}

// Root roots the tree by replacing trifurcation at the root by a
// bifurcation.
func (tree *Tree) Root(branchID int) error {
	if len(tree.Node.childNodes) == 2 {
		// no need to change
		return nil
	}

	reverse := false

	j := -1

	for i, node := range tree.Node.childNodes {
		if node.ID == branchID {
			j = i
			break
		}
	}
	if j < 0 {
		return errors.New("Rooting branch not found at the root node")
	}

	if j == 0 {
		reverse = true
	}

	newRoot := &Node{}
	if !reverse {
		newRoot.childNodes = append(newRoot.childNodes, tree.Node, tree.Node.childNodes[j])
	} else {
		newRoot.childNodes = append(newRoot.childNodes, tree.Node.childNodes[j], tree.Node)
	}

	tree.Node.BranchLength = tree.Node.childNodes[j].BranchLength / 2
	tree.Node.childNodes[j].BranchLength = tree.Node.BranchLength
	tree.Node.Class = tree.Node.childNodes[j].Class
	tree.Node.childNodes = append(tree.Node.childNodes[:j], tree.Node.childNodes[j+1:]...)
	tree.Node = newRoot

	newRoot.childNodes[0].Parent = newRoot
	newRoot.childNodes[1].Parent = newRoot

	// update cached parameters
	tree.nNodes++
	tree.nodeIDArray = nil
	tree.nodeOrder = nil

	return nil
}

// setLeafIDs labels all the terminal nodes of tree.
func (tree *Tree) setLeafIDs() {
	leafID := 0
	for node := range tree.Terminals() {
		node.LeafID = leafID
		leafID++
	}
}

// Node is a node of a tree.
type Node struct {
	// Name is the node label.
	Name string
	// BranchLength is the length of the  branch going to the node.
	BranchLength float64
	// Parent is the parent node (nil for the root).
	Parent *Node
	// ChildNodes is an array (slice) of children nodes.
	childNodes []*Node
	// Id is the node ID.
	ID int
	// LeafId is the leaf id.
	LeafID int
	// Class is a node class (specified by # symbol in newick).
	Class int
}

// NewNode creates a new Node.
func NewNode(parent *Node, nodeID int) (node *Node) {
	node = &Node{Parent: parent, ID: nodeID}
	return
}

// Copy creates copy of node with empty parent and children.
func (node *Node) Copy() *Node {
	return &Node{
		Name:         node.Name,
		BranchLength: node.BranchLength,
		childNodes:   make([]*Node, 0, len(node.childNodes)),
		ID:           node.ID,
		LeafID:       node.LeafID,
		Class:        node.Class,
	}
}

// AddChild adds a child node.
func (node *Node) AddChild(subNode *Node) {
	subNode.Parent = node
	node.childNodes = append(node.childNodes, subNode)
}

// BrString returns a newick string with branch labels for a node and
// its' subnodes.
func (node *Node) BrString() (s string) {
	if node.IsTerminal() {
		return fmt.Sprintf("%s*%d", node.Name, node.ID)
	}
	s += "("
	for i, child := range node.childNodes {
		s += child.BrString()
		if i != len(node.childNodes)-1 {
			s += ","
		}
	}
	s += fmt.Sprintf(")*%d", node.ID)
	if node.IsRoot() {
		s += ";"
	}
	return s
}

// ClassString returns a newick string with class labels for a node
// and its' subnodes.
func (node *Node) ClassString() (s string) {
	if node.IsTerminal() {
		if node.Class == 0 {
			return fmt.Sprintf("%s:%g", node.Name, node.BranchLength)
		}
		return fmt.Sprintf("%s#%d:%g", node.Name, node.Class, node.BranchLength)
	}
	s += "("
	for i, child := range node.childNodes {
		s += child.ClassString()
		if i != len(node.childNodes)-1 {
			s += ","
		}
	}
	s += ")"
	if node.Class != 0 {
		s += fmt.Sprintf("#%d", node.Class)
	}
	s += fmt.Sprintf(":%g", node.BranchLength)
	if node.IsRoot() {
		s += ";"
	}
	return s
}

// ShortClassString returns a newick string with class labels for a node
// and its' subnodes without branch lengths.
func (node *Node) ShortClassString() (s string) {
	if node.IsTerminal() {
		if node.Class == 0 {
			return fmt.Sprintf("%s", node.Name)
		}
		return fmt.Sprintf("%s#%d", node.Name, node.Class)
	}
	s += "("
	for i, child := range node.childNodes {
		s += child.ShortClassString()
		if i != len(node.childNodes)-1 {
			s += ","
		}
	}
	s += ")"
	if node.Class != 0 {
		s += fmt.Sprintf("#%d", node.Class)
	}
	if node.IsRoot() {
		s += ";"
	}
	return s
}

// String returns a newick string for a node and its' subnodes.
func (node *Node) String() (s string) {
	if node.IsTerminal() {
		return fmt.Sprintf("%s:%g", node.Name, node.BranchLength)
	}
	s += "("
	for i, child := range node.childNodes {
		s += child.String()
		if i != len(node.childNodes)-1 {
			s += ","
		}
	}
	s += fmt.Sprintf("):%g", node.BranchLength)
	if node.IsRoot() {
		s += ";"
	}
	return s
}

// LongString returns an extended string with all the ids, classes,
// etc for one node.
func (node *Node) LongString() (s string) {
	s = "<"
	if node.Parent == nil {
		s += "root, "
	}
	if node.Name != "" {
		s += "name=" + node.Name + ", "
	}
	s += fmt.Sprintf("Id=%v, BranchLength=%v", node.ID, node.BranchLength)
	if node.IsTerminal() {
		s += fmt.Sprintf(", TipId=%v", node.LeafID)
	}
	if node.Class != 0 {
		s += fmt.Sprintf(", Class=%v", node.Class)
	}
	s += ">"
	return
}

// FullString returns extended string for the whole tree.
func (node *Node) FullString() string {
	return strings.TrimSpace(node.prefixString(""))
}

// prefixString adds space prefix depending on the node level.
func (node *Node) prefixString(prefix string) (s string) {
	s = prefix + node.LongString() + "\n"
	for _, node := range node.childNodes {
		s += node.prefixString(prefix + "    ")
	}
	return
}

// ChildNodes returns a children nodes slice.
func (node *Node) ChildNodes() []*Node {
	return node.childNodes
}

// Walk traverses the tree nodes and sends them into a channel. An
// optional filtering function may be provided.
func (node *Node) Walk(ch chan *Node, filter func(*Node) bool) {
	if filter == nil || filter(node) {
		ch <- node
	}
	for _, node := range node.childNodes {
		node.Walk(ch, filter)
	}
}

// NSubNodes returns a total number of subnodes for a current node.
func (node *Node) NSubNodes() (size int) {
	for _, node := range node.childNodes {
		size += node.NSubNodes()
	}
	return size + 1
}

// IsRoot returns true if current node is a root.
func (node *Node) IsRoot() bool {
	return node.Parent == nil
}

// IsTerminal returns true if node is a leaf.
func (node *Node) IsTerminal() bool {
	return len(node.childNodes) == 0
}

// IsSpecial returns true if character is a special character
// (i.e. not identifier).
func IsSpecial(c rune) bool {
	switch c {
	case '(', ')', ':', '#', ';', ',':
		return true
	}
	return false

}

// NewickSplit is a function parsing a newick-formatted data.
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

// ParseNewick parses a newick formatted reader into the Tree type.
func ParseNewick(rd io.Reader) (tree *Tree, err error) {
	scanner := bufio.NewScanner(rd)

	scanner.Split(NewickSplit)

	nodeID := 0
	node := NewNode(nil, nodeID)
	tree = &Tree{Node: node}
	nodeID++

	mode := NORMAL

	for scanner.Scan() {
		text := scanner.Text()
		switch text {
		case "(":
			subNode := NewNode(nil, nodeID)
			nodeID++
			if node != nil {
				node.AddChild(subNode)
			}
			node = subNode

		case ",":
			if node.Parent == nil {
				return nil, errors.New("top level comma mismatch")
			}
			subNode := NewNode(nil, nodeID)
			nodeID++

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
			break
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
				node.Name = text
			}
		}
	}

	tree.setLeafIDs()

	return

}

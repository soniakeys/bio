package bio

import (
	"errors"
	"fmt"
	"strconv"
	"strings"
)

//--------------------------------------------------------------
// Stuff below the line will be replaced.  There are better ways.
// Prefix X added to mark things to be replaced.  Prefix changed to Z when a
// replacement is written.  When all X's are gone, Delete all the Z's

// XNewickNode is the node type used for both rooted and unrooted trees.
//
// Rooted trees are returned by ParseNewick and can represent any tree
// structure indicated by the Newick string representation.  A root node
// has no parent and its parent pointer must be nil.
// A rooted tree will have exactly one root node.
// The children of a node x must have their parent pointer set to x.
// A node with no children is a leaf.  There are no restrictions on
// the number of children a node may have.  There are no restrictions
// on which nodes must have names or weights.
//
// An unrooted tree is returned by XNewickNode.Unroot.  It is restricted
// to be a binary tree.  The parent pointer is not used to define tree
// structure.  Instead, structure is defined entirely with child pointers.
// An internal node will have exactly three children.
// A leaf node will will have exactly one child.  If node x has child y,
// then y must also have x as a child.
//
// If HasWeight is true and the node has a non-nil parent pointer, field
// "Weight" contains the weight or length of the arc from this node to
// its parent.  In the case of an unrooted tree, the parent pointer should
// be set appropriately to indicate which arc Weight applies to.  To weight
// all arcs in an unrooted tree, parent pointers will ultimately lead to
// a single internal node with no parent pointer.
//
// Field P is not part of the Newick format and can be used to associate
// arbitrary data with the node.
type XNewickNode struct {
	Name      string
	HasWeight bool
	Weight    float64
	Parent    *XNewickNode
	Children  []*XNewickNode
	P         interface{}
}

/*
func (n *XNewickNode) addChild(a *XNewickNode) {
	n.Children = append(n.Children, a)
}

func (n *XNewickNode) deleteChild(d *XNewickNode) (ok bool) {
	for i, ch := range n.Children {
		if ch == d {
			last := len(n.Children) - 1
			n.Children[i] = n.Children[last]
			n.Children = n.Children[:last]
			return true
		}
	}
	return false
}

func (n *XNewickNode) replaceChild(oldCh, newCh *XNewickNode) (ok bool) {
	for i, ch := range n.Children {
		if ch == oldCh {
			n.Children[i] = newCh
			return true
		}
	}
	return false
}
*/

// XNewickNames represents an index from node names to node pointers.
type XNewickNames map[string]*XNewickNode

// RootedDistance returns the distance between two named nodes in a rooted
// tree.
//
// The reciever must be an index into a rooted tree.
// The distance returned is the number of edges separating the named nodes.
func (names XNewickNames) RootedDistance(name1, name2 string) (int, error) {
	n1, ok := names[name1]
	if !ok {
		return 0, fmt.Errorf("%s not in tree", name1)
	}
	n2, ok := names[name2]
	if !ok {
		return 0, fmt.Errorf("%s not in tree", name2)
	}
	if name1 == name2 {
		return 0, nil
	}
	walked := map[*XNewickNode]int{}
	walked[n1] = 0
	walked[n2] = 0
	for d := 1; ; d++ {
		if n1.Parent != nil {
			if d2, ok := walked[n1.Parent]; ok {
				// parent is on path 2
				return d + d2, nil
			}
			n1 = n1.Parent
			walked[n1] = d
		}
		if n2.Parent != nil {
			if d1, ok := walked[n2.Parent]; ok {
				// parent is on path1
				return d + d1, nil
			}
			n2 = n2.Parent
			walked[n2] = d
		}
	}
	return 0, fmt.Errorf("invalid tree")
}

// RootedDistanceByWeight returns the distance between two named nodes in
// a rooted tree.
//
// The reciever must be an index into a rooted tree with edge weights
// The distance returned is the sum of the weights of the edge separating
// the named nodes.
func (names XNewickNames) RootedDistanceByWeight(name1, name2 string) (float64, error) {
	n1, ok := names[name1]
	if !ok {
		return 0, fmt.Errorf("%s not in tree", name1)
	}
	n2, ok := names[name2]
	if !ok {
		return 0, fmt.Errorf("%s not in tree", name2)
	}
	if name1 == name2 {
		return 0, nil
	}
	part1 := map[*XNewickNode]float64{n1: 0}
	part2 := map[*XNewickNode]float64{n2: 0}
	var all1, all2 float64
	for {
		if n1.Parent != nil {
			if !n1.HasWeight {
				name1 := n1.Name
				if name1 == "" {
					name1 = "(unnamed)"
				}
				name2 := n1.Parent.Name
				if name2 == "" {
					name2 = "(unnamed)"
				}
				return 0, fmt.Errorf("arc %s to %s is unweighted", name1, name2)
			}
			all1 += n1.Weight
			n1 = n1.Parent
			if p2, ok := part2[n1]; ok {
				// n1 is on path 2
				return all1 + p2, nil
			}
			part1[n1] = all1
		}
		if n2.Parent != nil {
			if !n2.HasWeight {
				name1 := n2.Name
				if name1 == "" {
					name1 = "(unnamed)"
				}
				name2 := n2.Parent.Name
				if name2 == "" {
					name2 = "(unnamed)"
				}
				return 0, fmt.Errorf("arc %s to %s is unweighted", name1, name2)
			}
			all2 += n2.Weight
			n2 = n2.Parent
			if p1, ok := part1[n2]; ok {
				// n2 is on path 1
				return all2 + p1, nil
			}
			part2[n2] = all2
		}
	}
	return 0, fmt.Errorf("invalid tree")
}

// ParseNewick parses the string representation of Newick data and returns
// the root node of a rooted tree.
//
// Argument s must have a terminanting semicolon.  There can be nothing but
// whitespace follwing the semicolon.
// If argument "names" is non-nil, it is populated as well.
func XParseNewick(s string, names XNewickNames) (*XNewickNode, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return nil, errors.New("no data")
	}
	last := len(s) - 1
	if s[last] != ';' {
		return nil, errors.New("string not terminated with ;")
	}
	if len(s) == 1 {
		return nil, nil // empty tree is valid
	}
	if names != nil {
		for k := range names { // clear the map
			delete(names, k)
		}
	}
	p := &parser{rem: s[:last], names: names}
	p.gettok()
	n, err := p.parseSubtree()
	if err != nil {
		return nil, err
	}
	if p.rem > "" {
		fmt.Println("debug: p.rem:", p.rem)
		return nil, errors.New("unparsed text follows complete tree")
	}
	return n, nil
}

type parser struct {
	rem   string
	tok   string
	names XNewickNames
}

func (p *parser) gettok() {
	if p.rem == "" {
		p.tok = ""
		return
	}
	switch p.rem[0] {
	case '(', ')', ',':
		p.tok = string(p.rem[0])
		p.rem = strings.TrimSpace(p.rem[1:])
		return
	}
	if x := strings.IndexAny(p.rem, "(),"); x > 0 {
		p.tok = strings.TrimSpace(p.rem[:x])
		p.rem = p.rem[x:]
	} else {
		p.tok = p.rem
		p.rem = ""
	}
}

func (p *parser) parseSubtree() (n *XNewickNode, err error) {
	if p.tok == "(" {
		// internal node
		return p.parseSet()
	}
	// leaf node
	n = &XNewickNode{}
	if p.tok != ")" && p.tok != "," {
		err = p.nameWeight(n)
	}
	return
}

func (p *parser) parseSet() (*XNewickNode, error) {
	n := &XNewickNode{}
	p.gettok() // get token after (
	for {
		child, err := p.parseSubtree()
		if err != nil {
			return n, err
		}
		n.Children = append(n.Children, child)
		child.Parent = n
		if p.tok != "," {
			break
		}
		p.gettok()
	}
	if p.tok != ")" {
		return nil, errors.New("expected )")
	}
	p.gettok()
	switch p.tok {
	case ")", ",", "(":
		return n, nil
	}
	if err := p.nameWeight(n); err != nil {
		return nil, err
	}
	return n, nil
}

// add name and weight to either internal or leaf node
func (p *parser) nameWeight(n *XNewickNode) (err error) {
	w := strings.Index(p.tok, ":")
	if w < 0 {
		n.Name = p.tok
	} else {
		if n.Weight, err = strconv.ParseFloat(p.tok[w+1:], 64); err != nil {
			return err
		}
		n.HasWeight = true
		n.Name = p.tok[:w]
	}
	if p.names != nil && n.Name > "" {
		if _, exists := p.names[n.Name]; exists {
			return fmt.Errorf("name %s appears more than once", n.Name)
		}
		p.names[n.Name] = n
	}
	p.gettok() // get token after name:weight
	return nil
}

// String satisfies fmt.Stringer by formatting a tree of XNewickNodes into
// the Newick text format.
//
// Only child nodes are followed.  Passing the root of a rooted tree will
// format the entire tree.  Passing a node under the root will format only
// that sub tree.  In the case of an unrooted tree, any node may be passed
// and String will format the entire tree rooted at than node.
func (n *XNewickNode) String() string {
	if n == nil {
		return ";"
	}
	return rNewickNodeString(nil, n) + ";"
}

// recursive stringer
func rNewickNodeString(from, to *XNewickNode) string {
	s := ""
	if len(to.Children) > 0 {
		s += "("
		one := false
		for _, ch := range to.Children {
			if ch != from {
				if one {
					s += ","
				}
				s += rNewickNodeString(to, ch)
				one = true
			}
		}
		s += ")"
	}
	s += to.Name
	if to.HasWeight {
		s += fmt.Sprintf(":%g", to.Weight)
	}
	return s
}

// Unroot traverses rooted input and constructs unrooted output.
//
// Rooted input has the shape returned by ParseNewick, which retains
// the shape of the input Newick syntax.  It not only has a node identified
// as root, but the links between nodes are asymetric; there is a parent/
// child relationship.
//
// Unrooted output has no single node identified as root.  Instead it
// returns an index to the leaves which can be used to enter the tree
// from any leaf.  Further, the only parent/child asymetry is between
// internal nodes and leaf nodes.  Among the internal nodes, links are
// symetric and use a child/child relationship.
//
// Unrooted output is validated to be a binary tree.
// Internal nodes must be unnamed, have 3 children and no parent.
// Leaf nodes must be named, have a parent and no children.
// The names must be unique and are returned as a XNewickNames map.
// The map contains pointers to leaf nodes and the tree can be traversed
// from there.
//
// Incomplete:  Not all cases are handled.  Too restrictive.  More
// cases can be represented.
func (root *XNewickNode) Unrooted() (XNewickNames, error) {
	names := XNewickNames{}
	// the case of the CTBL test data format:  a named root with one child.
	// create node that will become leaf, convert tree as a child of this node,
	// then move child to parent field.
	if root.Name > "" && len(root.Children) == 1 {
		uLeaf := &XNewickNode{Name: root.Name}
		names[root.Name] = uLeaf
		if err := rUnrooted(root.Children[0], uLeaf, names); err != nil {
			return nil, err
		}
		uLeaf.Parent = uLeaf.Children[0]
		uLeaf.Children = nil
		return names, nil
	}
	// otherwise
	return nil, errors.New("unhandled input format")
}

// recursive parser.
// d for the input directed node, u for result undirected node
func rUnrooted(d *XNewickNode, uIn *XNewickNode, names XNewickNames) error {
	if len(d.Children) == 0 {
		// leaf node
		if d.Name == "" {
			return errors.New("unnamed leaf")
		}
		if _, exists := names[d.Name]; exists {
			return fmt.Errorf("name %s not unique", d.Name)
		}
		uOut := &XNewickNode{Name: d.Name, Parent: uIn}
		uIn.Children = append(uIn.Children, uOut)
		names[d.Name] = uOut
		return nil
	}
	// internal node
	if len(d.Children) != 2 {
		return fmt.Errorf("internal node with %d children", len(d.Children))
	}
	if d.Name > "" {
		return fmt.Errorf("named internal node: %s", d.Name)
	}
	uOut := &XNewickNode{Children: []*XNewickNode{uIn}}
	uIn.Children = append(uIn.Children, uOut)
	if err := rUnrooted(d.Children[0], uOut, names); err != nil {
		return err
	}
	if err := rUnrooted(d.Children[1], uOut, names); err != nil {
		return err
	}
	return nil
}

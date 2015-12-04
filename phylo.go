package bio

import (
	"errors"
	"fmt"
	"strconv"
	"strings"

	"github.com/soniakeys/graph"
)

// PhyloList represents a rooted tree using a graph.FromList.
//
// It has node names and edge weights as in the Newick format.
//
// This is a compact minimimal representation but may not be convenient
// for some algorithms.  See PhyloRootedTree for an alternative representation.
type PhyloList struct {
	Tree       graph.FromList // tree structure
	Nodes      []PhyloNode    // parallel to Tree
	Root       int            // root node of tree
	NumNames   int            // count of Name > ""
	NumWeights int            // count of HasWt == true
}

// PhyloRootedTree represents a rooted tree.
//
// It has node names and edge weights as in the Newick format.
type PhyloRootedTree struct {
	Tree       graph.LabeledAdjacencyList // tree structure
	Nodes      []PhyloNode                // parallel to Tree
	Root       int                        // root node of tree
	NumNames   int                        // count of Name > ""
	NumWeights int                        // count of HasWt == true
}

// PhyloNode represents data for a single node of a rooted tree.
type PhyloNode struct {
	Name      string
	HasWeight bool
	Weight    float64 // weight of edge to parent node.
}

func (t *PhyloList) RootedTree() *PhyloRootedTree {
	return &PhyloRootedTree{
		Tree:       t.Tree.TransposeLabeled(),
		Nodes:      t.Nodes,
		Root:       t.Root,
		NumNames:   t.NumNames,
		NumWeights: t.NumWeights,
	}
}

// Newick serialializes to Newick format.
func (t *PhyloRootedTree) Newick() string {
	tr := t.Tree
	var f func(graph.Half) string
	f = func(p graph.Half) (s string) {
		to := tr[p.To]
		if len(to) > 0 { // format children
			c := make([]string, len(to))
			for i, to := range to {
				c[i] = f(to)
			}
			s = fmt.Sprintf("(%s)", strings.Join(c, ","))
		}
		s += t.Nodes[p.To].Name
		if p.Label >= 0 {
			if nd := t.Nodes[p.Label]; nd.HasWeight {
				s = fmt.Sprintf("%s:%g", s, nd.Weight)
			}
		}
		return s
	}
	return f(graph.Half{To: t.Root, Label: -1}) + ";"
}

type newickParser struct {
	rem string
	tok string
	rt  *PhyloList
}

// ParseNewick parses Newick format.
//
// Argument s must have a terminating semicolon.  There can be nothing but
// whitespace follwing the semicolon.
func ParseNewick(s string) (*PhyloList, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return nil, errors.New("no data")
	}
	last := len(s) - 1
	if s[last] != ';' {
		return nil, errors.New("string not terminated with ;")
	}
	rt := &PhyloList{} // zero value is valid empty tree
	if len(s) == 1 {
		return rt, nil // empty tree
	}
	// tree is not empty, create root
	rt.Tree.Paths = []graph.PathEnd{{From: -1, Len: 1}}
	rt.Tree.MaxLen = 1
	rt.Nodes = []PhyloNode{{}}

	p := &newickParser{rem: s[:last], rt: rt}
	p.gettok()
	if err := p.parseSubtree(0); err != nil {
		return nil, err
	}
	if p.rem > "" {
		if len(p.rem) > 30 {
			p.rem = p.rem[:27] + "..."
		}
		// fmt.Println("debug: p.rem:", p.rem)
		return nil, errors.New("unparsed text follows complete tree: " + p.rem)
	}
	return rt, nil
}

func (p *newickParser) gettok() {
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

func (p *newickParser) parseSubtree(n int) (err error) {
	if p.tok == "(" {
		// internal node
		return p.parseSet(n)
	}
	// leaf node
	p.rt.Tree.Leaves.SetBit(&p.rt.Tree.Leaves, n, 1)
	if p.tok != ")" && p.tok != "," {
		err = p.nameWeight(n)
	}
	return
}

// add name and weight to node
func (p *newickParser) nameWeight(n int) (err error) {
	pn := &p.rt.Nodes[n]
	tok := p.tok
	if w := strings.Index(tok, ":"); w >= 0 {
		if pn.Weight, err = strconv.ParseFloat(tok[w+1:], 64); err != nil {
			return err
		}
		pn.HasWeight = true
		p.rt.NumWeights++
		tok = tok[:w]
	}
	if tok > "" {
		pn.Name = tok
		p.rt.NumNames++
	}
	p.gettok() // get token after name:weight
	return nil
}

func (p *newickParser) parseSet(n int) error {
	p.gettok()                  // get token after (
	fl := &p.rt.Tree            // dereference FromList
	pLen := fl.Paths[n].Len + 1 // path length to nodes in this set
	if pLen > fl.MaxLen {
		fl.MaxLen = pLen
	}
	for {
		// create child node
		cn := len(fl.Paths)
		fl.Paths = append(fl.Paths, graph.PathEnd{From: n, Len: pLen})
		p.rt.Nodes = append(p.rt.Nodes, PhyloNode{})

		if err := p.parseSubtree(cn); err != nil {
			return err
		}
		if p.tok != "," {
			break
		}
		p.gettok()
	}
	if p.tok != ")" {
		return errors.New("expected )")
	}
	p.gettok()
	switch p.tok {
	case ")", ",", "(":
		return nil
	}
	if err := p.nameWeight(n); err != nil {
		return err
	}
	return nil
}

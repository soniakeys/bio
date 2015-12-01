package bio

import (
	"errors"
	"fmt"
	"strconv"
	"strings"

	"github.com/soniakeys/graph"
)

// NewickTree represents data parsed from Newick format.
type NewickTree struct {
	Tree       graph.FromList // tree structure
	Nodes      []NewickNode   // parallel to Tree
	NumNames   int            // count of Name > ""
	NumWeights int            // count of HasWt == true
}

// NewickNode represents data for a single node of a Newick tree.
type NewickNode struct {
	Name      string
	HasWeight bool
	Weight    float64
}

// String serialializes to Newick format.
func (t *NewickTree) String() string {
	m := map[int]string{}
	for n := t.Tree.Leaves.BitLen(); n > 0; {
		n--
		if t.Tree.Leaves.Bit(n) == 1 {
			nd := &t.Nodes[n]
			fr := t.Tree.Paths[n].From
			s := fmt.Sprintf("%s,%s", m[fr], nd.Name)
			if nd.HasWeight {
				s = fmt.Sprintf("%s:%g", s, nd.Weight)
			}
			m[fr] = s
		}
	}
	for {
		for n, c := range m {
			if n == -1 {
				return c[1:] + ";"
			}
			nd := &t.Nodes[n]
			fr := t.Tree.Paths[n].From
			s := fmt.Sprintf("%s,(%s)%s", m[fr], c[1:], nd.Name)
			if nd.HasWeight {
				s = fmt.Sprintf("%s:%g", s, nd.Weight)
			}
			m[fr] = s
		}
	}
}

type newickParser struct {
	rem string
	tok string
	nt  *NewickTree
}

// ParseNewick parses Newick format.
//
// Argument s must have a terminating semicolon.  There can be nothing but
// whitespace follwing the semicolon.
func ParseNewick(s string) (*NewickTree, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return nil, errors.New("no data")
	}
	last := len(s) - 1
	if s[last] != ';' {
		return nil, errors.New("string not terminated with ;")
	}
	nt := &NewickTree{} // zero value is valid empty tree
	if len(s) == 1 {
		return nt, nil // empty tree
	}
	// tree is not empty, create root
	nt.Tree.Paths = []graph.PathEnd{{From: -1, Len: 1}}
	nt.Tree.MaxLen = 1
	nt.Nodes = []NewickNode{{}}

	p := &newickParser{rem: s[:last], nt: nt}
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
	return nt, nil
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
	p.nt.Tree.Leaves.SetBit(&p.nt.Tree.Leaves, n, 1)
	if p.tok != ")" && p.tok != "," {
		err = p.nameWeight(n)
	}
	return
}

// add name and weight to node
func (p *newickParser) nameWeight(n int) (err error) {
	pn := &p.nt.Nodes[n]
	tok := p.tok
	if w := strings.Index(tok, ":"); w >= 0 {
		if pn.Weight, err = strconv.ParseFloat(tok[w+1:], 64); err != nil {
			return err
		}
		pn.HasWeight = true
		p.nt.NumWeights++
		tok = tok[:w]
	}
	if tok > "" {
		pn.Name = tok
		p.nt.NumNames++
	}
	p.gettok() // get token after name:weight
	return nil
}

func (p *newickParser) parseSet(n int) error {
	p.gettok()                  // get token after (
	fl := &p.nt.Tree            // dereference FromList
	pLen := fl.Paths[n].Len + 1 // path length to nodes in this set
	if pLen > fl.MaxLen {
		fl.MaxLen = pLen
	}
	for {
		// create child node
		cn := len(fl.Paths)
		fl.Paths = append(fl.Paths, graph.PathEnd{From: n, Len: pLen})
		p.nt.Nodes = append(p.nt.Nodes, NewickNode{})

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

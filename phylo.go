package bio

import (
	"errors"
	"fmt"
	"math"
	"math/big"
	"strconv"
	"strings"

	"github.com/soniakeys/graph"
)

// PhyloList represents a rooted tree using a graph.FromList.
//
// It has node names and arc weights as in the Newick format.
//
// This is a compact and fairly minimimal representation that still
// facilites computation.  In this representation arcs are directed
// from the leaves toward the root.
//
// Some algorithms require a different representation however.
// See PhyloRootedTree for an alternative representation with a different
// method set.
type PhyloList struct {
	List       graph.FromList    // tree encoded as a parent list
	Root       graph.NI          // root node of tree
	Nodes      []PhyloRootedNode // parallel to Tree
	NumLeaves  int               // count of leaf nodes
	NumNames   int               // count of Name > ""
	NumWeights int               // count of HasWt == true
}

// PhyloRootedTree represents a rooted tree.
//
// It has node names and arc weights as in the Newick format.
//
// In this representation arcs are directed from the root toward the leaves.
type PhyloRootedTree struct {
	Tree       graph.Directed    // adjacency list tree structure
	Root       graph.NI          // root node of tree
	Nodes      []PhyloRootedNode // parallel to Tree
	NumLeaves  int               // count of leaf nodes
	NumNames   int               // count of Name > ""
	NumWeights int               // count of HasWt == true
}

// PhyloRootedNode represents data for a single node of a rooted tree.
type PhyloRootedNode struct {
	Name      string
	HasWeight bool
	Weight    float64 // weight of arc from parent node.
}

// RootedTree construsts a PhyloRootedTree equivalent to a PhyloList.
func (l *PhyloList) RootedTree() *PhyloRootedTree {
	return &PhyloRootedTree{
		Tree:       l.List.Transpose(),
		Nodes:      l.Nodes,
		Root:       l.Root,
		NumLeaves:  l.NumLeaves,
		NumNames:   l.NumNames,
		NumWeights: l.NumWeights,
	}
}

// List constructs a PhyloList equivalent to a PhyloRootedTree.
func (t *PhyloRootedTree) List() *PhyloList {
	l, n := t.Tree.FromList()
	if n < 0 {
		return nil
	}
	return &PhyloList{
		List:       *l,
		Nodes:      t.Nodes,
		Root:       t.Root,
		NumLeaves:  t.NumLeaves,
		NumNames:   t.NumNames,
		NumWeights: t.NumWeights,
	}
}

// Newick serializes to Newick format.
func (t *PhyloRootedTree) Newick() string {
	tr := t.Tree
	var f func(graph.NI) string
	f = func(p graph.NI) (s string) {
		to := tr.AdjacencyList[p]
		if len(to) > 0 { // format children
			c := make([]string, len(to))
			for i, to := range to {
				c[i] = f(to)
			}
			s = fmt.Sprintf("(%s)", strings.Join(c, ","))
		}
		nd := t.Nodes[p]
		s += nd.Name
		if nd.HasWeight {
			s = fmt.Sprintf("%s:%g", s, nd.Weight)
		}
		return s
	}
	return f(t.Root) + ";"
}

type newickParser struct {
	rem string
	tok string
	pl  *PhyloList
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
	pl := &PhyloList{} // zero value is valid empty tree
	if len(s) == 1 {
		return pl, nil // empty tree
	}
	// tree is not empty, create root
	pl.List.Paths = []graph.PathEnd{{From: -1, Len: 1}}
	pl.List.MaxLen = 1
	pl.Nodes = []PhyloRootedNode{{}}

	p := &newickParser{rem: s[:last], pl: pl}
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
	return pl, nil
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

func (p *newickParser) parseSubtree(n graph.NI) (err error) {
	if p.tok == "(" {
		// internal node
		return p.parseSet(n)
	}
	// leaf node
	p.pl.NumLeaves++
	p.pl.List.Leaves.SetBit(&p.pl.List.Leaves, int(n), 1)
	if p.tok != ")" && p.tok != "," {
		err = p.nameWeight(n)
	}
	return
}

// add name and weight to node
func (p *newickParser) nameWeight(n graph.NI) (err error) {
	pn := &p.pl.Nodes[n]
	tok := p.tok
	if w := strings.Index(tok, ":"); w >= 0 {
		if pn.Weight, err = strconv.ParseFloat(tok[w+1:], 64); err != nil {
			return err
		}
		pn.HasWeight = true
		p.pl.NumWeights++
		tok = tok[:w]
	}
	if tok > "" {
		pn.Name = tok
		p.pl.NumNames++
	}
	p.gettok() // get token after name:weight
	return nil
}

func (p *newickParser) parseSet(n graph.NI) error {
	p.gettok()                  // get token after (
	fl := &p.pl.List            // dereference FromList
	pLen := fl.Paths[n].Len + 1 // path length to nodes in this set
	if pLen > fl.MaxLen {
		fl.MaxLen = pLen
	}
	for {
		// create child node
		cn := graph.NI(len(fl.Paths))
		fl.Paths = append(fl.Paths, graph.PathEnd{From: n, Len: pLen})
		p.pl.Nodes = append(p.pl.Nodes, PhyloRootedNode{})

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

// PathLen returns the number of arcs between the nodes.
func (l *PhyloList) PathLen(a, b graph.NI) int {
	p := l.List.Paths
	return p[a].Len + p[b].Len - 2*p[l.List.CommonAncestor(a, b)].Len
}

// Distance returns the sum of arc weights between two nodes
func (l *PhyloList) Distance(a, b graph.NI) (d float64) {
	// code similar to graph.CommonAncestor
	p := l.List.Paths
	n := l.Nodes
	if a < 0 || b < 0 || a >= graph.NI(len(p)) || b >= graph.NI(len(p)) {
		return math.NaN()
	}
	if p[a].Len < p[b].Len {
		a, b = b, a
	}
	for bl := p[b].Len; p[a].Len > bl; {
		d += n[a].Weight
		a = p[a].From
	}
	for a != b {
		d += n[a].Weight + n[b].Weight
		a = p[a].From
		b = p[b].From
	}
	return d
}

// NodeMap constructs a map from node names to node indexes.
//
// Only named nodes are included in the map and names are assumed to be
// unique.  That is, each distinct non-empty name will map to a single node
// index.
func (l *PhyloList) NodeMap() map[string]graph.NI {
	m := map[string]graph.NI{}
	for n, nd := range l.Nodes {
		if nd.Name > "" {
			m[nd.Name] = graph.NI(n)
		}
	}
	return m
}

// CharacterTable lists non-trivial characters of a phylogenetic tree.
//
// A character is represented as bits of a big.Int with bits set by node number
// of the tree.  A non-trivial character corresponds to an internal arc of
// the tree; that is, an arc not connecting a leaf node.  For a single
// character, bits of the big.Int will be set to 1 for nodes of the subtree
// pointed to by the arc, other bits will be 0.
func (t *PhyloRootedTree) CharacterTable() []big.Int {
	g := t.Tree
	chars := make([]big.Int, 0, len(g.AdjacencyList)-3)
	var f func(graph.NI) big.Int
	f = func(p graph.NI) big.Int {
		var pb big.Int
		pb.SetBit(&pb, int(p), 1)
		for _, to := range g.AdjacencyList[p] {
			tb := f(to)
			if len(g.AdjacencyList[to]) > 0 {
				chars = append(chars, tb)
			}
			pb.Or(&pb, &tb)
		}
		return pb
	}
	n := t.Root
	if to := g.AdjacencyList[n]; len(to) == 1 {
		n = to[0] // skip root leaf
	}
	f(n)
	return chars
}

// CharacterTable produces a phylogenetic character table from
// a set of equal length symbol strings.
//
// A character is considered to exist for symbol position when at least
// two strings have the the same symbol and at least two two strings have
// a different symbol.  In this case a modal character is determined and
// strings with a symbol other than the nominal modal will have the bit
// set to 1 in the corresponding position of the character.
//
// At least four strings are required.  An error is returned for < 4 strings
// or for strings of unequal length.
func (sk StrKmers) CharacterTable() (ct []big.Int, pos []int, err error) {
	if len(sk) < 4 {
		return nil, nil, errors.New("not enough strings to characterize")
	}
	// validate string lengths equal
	s0 := sk[0]
	sLen := len(sk[0])
	if sLen == 0 {
		return nil, nil, errors.New("can't characterize empty strings")
	}
	for _, s := range sk[1:] {
		if len(s) != sLen {
			return nil, nil, errors.New("strings different lengths")
		}
	}
	for i := range s0 {
		// loop1, find mode
		var f [256]int
		var mode int
		var modal byte
		for _, s := range sk {
			sym := s[i]
			f[sym]++
			if f[sym] > mode {
				mode = f[sym]
				modal = sym
			}
		}
		if len(sk)-mode < 2 {
			continue // no or trivial split
		}
		// loop 2, add character
		var c big.Int
		for p, s := range sk {
			if s[i] != modal {
				c.SetBit(&c, p, 1)
			}
		}
		ct = append(ct, c)
		pos = append(pos, i)
	}
	return
}

// MaxParsimony implements "small parsimony" -- max parsimony for
// the given tree structure.
//
// MaxParsimony computes kmers for internal nodes of the tree.  It solves
// for kmers that minimize the total tree weight, as the sum of all arc
// weights, where an arc weight is the DNA8 hamming distance between a node
// and its parent.
//
// The argument kmers is parallel to the nodes of t.  On entry, all kmers
// must be allocated and kmers for leaf nodes must have valid DNA8 content.
//
// MaxParsimonly fills in the content for kmers of internal nodes, and assigns
// weights to all nodes of the tree.  The root will be assigned weight 0,
// but with HasWeight set to false.
func (t *PhyloRootedTree) MaxParsimony(kmers Kmers) (total int) {
	tree := t.Tree.AdjacencyList
	nodes := t.Nodes
	cost4 := make([][4]float64, len(tree)) // (dis)parsimony scores
	inf4 := [4]float64{math.Inf(1), math.Inf(1), math.Inf(1), math.Inf(1)}

	var score func(graph.NI, int)
	score = func(n graph.NI, sx int) {
		if len(tree[n]) == 0 {
			copy(cost4[n][:], inf4[:])
			cost4[n][kmers[n][sx]>>1&3] = 0
			return
		}
		//fmt.Println("scoring node", n)
		for _, ch := range tree[n] {
			score(ch, sx)
		}
		for k := range "ACTG" {
			//fmt.Println("  base index", k)
			d := 0. // s(k)(v) from the text, the (dis)pars score for sym k
			for _, ch := range tree[n] {
				//fmt.Println("    child", ch, "dis[ch]", dis[ch])
				min := math.Inf(1)
				for i, dch := range cost4[ch] {
					if i != k {
						dch++
					}
					if dch < min {
						min = dch
					}
				}
				d += min
			}
			cost4[n][k] = d
		}
		//fmt.Println("node", n, "scored:", dis[n])
	}

	var labelNodes func(graph.NI, int, int)
	labelNodes = func(n graph.NI, sx, axp int) {
		if len(tree[n]) == 0 { // leaf nodes come with labels
			// just need to note distance
			if int(kmers[n][sx]>>1&3) != axp {
				nodes[n].Weight++
				total++
			}
			return
		}
		min := math.Inf(1)
		var xMin int
		for x, d := range cost4[n] {
			if x != axp {
				d++
			}
			if d < min {
				min = d
				xMin = x
			}
		}
		kmers[n][sx] = "ACTG"[xMin]
		if xMin != axp {
			nodes[n].Weight++
			total++
		}
		for _, ch := range tree[n] {
			labelNodes(ch, sx, xMin)
		}
	}

	// initialize
	for n := range nodes {
		nodes[n].Weight = 0
		nodes[n].HasWeight = true
	}
	// outer loop goes by symbol position
	for sx := range kmers[0] {
		// depth first pass 1: accumulate dis-parsimony scores
		score(t.Root, sx)
		// depth first pass 2: build labels for internal nodes
		labelNodes(t.Root, sx, -1)
	}
	total -= len(kmers[0])
	nodes[t.Root].Weight = 0
	nodes[t.Root].HasWeight = false
	t.NumWeights = len(tree) - 1
	return
}

// DistanceMatrix constructs a distance matrix from a PhyloList.
func (l *PhyloList) DistanceMatrix(leaves []graph.NI) [][]float64 {
	m := make([][]float64, len(leaves))
	for i := range m {
		m[i] = make([]float64, len(leaves))
	}
	for i := 1; i < len(leaves); i++ {
		li := leaves[i]
		mi := m[i]
		for j, lj := range leaves[:i] {
			d := l.Distance(li, lj)
			mi[j] = d
			m[j][i] = d
		}
	}
	return m
}

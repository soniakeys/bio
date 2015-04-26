package bio

import (
	"math"
	"math/big"

	"github.com/soniakeys/graph"
)

// utree.go
//
// Functions for unrooted binary trees, including maximum parsimony.
// Currently these are all trees of DNA8 kmers.  They use types from
// the graph library.  One function takes a tree as a directed graph,
// all others take trees as edge-labeled undirected graphs.

// UTree is an unrooted binary tree with DNA8 labeled nodes, and float64
// cost labeled edges.
type UTree struct {
	NLeaves int
	T       graph.LabeledAdjacencyList // len(T) is 2*NLeaves-2
	Kmers   Kmers                      // len(Kmers) is len(T)
	Costs   []float64                  // edge costs, len(Costs) = len(T)-1
}

// but first here's the rooted tree function,

// DNA8MaxParsimonyRooted solves kmers and costs for maximum parsimony
// (minimum total cost) in a rooted phylogenic tree.
//
// The argument tree is a directed graph, nodes 0:len(leaves) must be leaf
// nodes, node len(tree)-1 must be the root.
//
// Returned kmers will be the leaves argument with appended kmers labeling
// internal nodes and the root.  Returned cost slice is parallel to kmers
// (and the input tree).  Cost for a node is the hamming distance from the
// node's parent.  Cost 0 is returned for the root.
func DNA8MaxParsimonyRooted(tree graph.AdjacencyList, leaves Kmers) (kmers Kmers, costs []float64) {
	root := len(tree) - 1
	m := len(leaves[0])
	kmers = append(leaves, make(Kmers, len(tree)-len(leaves))...)
	for i := len(leaves); i < len(tree); i++ {
		kmers[i] = make(DNA8, m)
	}
	costs = make([]float64, len(tree))     // cost is edge length from parent
	cost4 := make([][4]float64, len(tree)) // (dis)parsimony scores
	inf4 := [4]float64{math.Inf(1), math.Inf(1), math.Inf(1), math.Inf(1)}
	var score func(int, int)
	score = func(n, sx int) {
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

	var labelNodes func(int, int, int)
	labelNodes = func(n, sx, axp int) {
		if len(tree[n]) == 0 { // leaf nodes come with labels
			// just need to note distance
			if int(kmers[n][sx]>>1&3) != axp {
				costs[n]++
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
			costs[n]++
		}
		for _, ch := range tree[n] {
			labelNodes(ch, sx, xMin)
		}
	}

	// outer loop goes by symbol position
	for sx := range kmers[0] {
		// depth first pass 1: accumulate dis-parsimony scores
		score(root, sx)
		// depth first pass 2: build labels for internal nodes
		labelNodes(root, sx, -1)
	}
	costs[root] = 0
	return
}

// dna8MaxParsimonyUnrooted uses above DNA8MaxParsimonyRooted as a subroutine.
// Thus it builds directed trees from the input undirected tree.
// The corresponding exported method on UTree should be equivalent but more
// efficient since it works directly on the undirected tree.  This version
// left for reference and testing.
//
// dna8MaxParsimonyUnrooted solves kmers and costs for maximum parsimony
// (minimum total cost) in a rooted phylogenic tree.
//
// The argument tree is a undirected graph, nodes 0:len(leaves) must be leaf
// nodes.  It must be a tree and so have len(tree)-1 edges, labeled with
// numbers 0:len(tree)-1.
//
// Returned kmers will be the leaves argument with appended kmers labeling
// internal nodes and the root.  Returned cost slice is indexed by the
// graph edge labels.
func dna8MaxParsimonyUnrooted(tree graph.LabeledAdjacencyList, leaves Kmers) (kmers Kmers, costs []float64) {
	// for each edge, construct rooted tree, call DNA8MaxParsimonyRooted().
	// df traversal to get each edge just once.
	root := len(tree) // new node
	rt := make(graph.AdjacencyList, root+1)
	// store edges labels to with associate costs in the end.
	labels := make([]int, len(tree))
	var pop big.Int
	var populate func(int, graph.Half)
	populate = func(parent int, ch graph.Half) {
		pop.SetBit(&pop, parent, 1)
		rt[parent] = append(rt[parent], ch.To)
		labels[ch.To] = ch.Label
		parent = ch.To
		for _, ch := range tree[parent] {
			if pop.Bit(ch.To) == 0 {
				populate(parent, ch)
			}
		}
	}
	split := func(n1 int, n2 graph.Half) {
		pop = big.Int{}
		pop.SetBit(&pop, n1, 1)
		pop.SetBit(&pop, n2.To, 1)
		populate(root, graph.Half{n1, n2.Label})
		populate(root, n2)
	}

	var minC = math.Inf(1)
	var minTree = make([][]int, len(rt))
	var minKmers Kmers
	var minCosts []float64
	var minEdgeLabels = make([]int, len(labels))

	var vis big.Int
	var rootEdges func(int)
	rootEdges = func(n int) {
		vis.SetBit(&vis, n, 1)
		for _, ch := range tree[n] {
			if vis.Bit(ch.To) == 1 {
				continue
			}
			// root tree at this edge.
			for i := range rt {
				rt[i] = nil // clear tree
			}
			split(n, ch)
			kmers, costs := DNA8MaxParsimonyRooted(rt, leaves)
			c := 0.
			for _, ec := range costs {
				c += ec
			}
			if c < minC {
				minC = c
				copy(minTree, rt)
				minKmers = kmers
				minCosts = costs
				copy(minEdgeLabels, labels)
			}
			rootEdges(ch.To) // df traversal to all edges
		}
	}
	rootEdges(0)
	costs = make([]float64, len(tree)-1)
	for n, l := range minEdgeLabels {
		costs[l] += minCosts[n]
	}
	return minKmers[:len(tree)], costs
}

func (u *UTree) MaxParsimonyEdge(from int, to graph.Half) float64 {
	cost4 := make([][4]float64, len(u.T)) // (dis)parsimony scores
	//	var vis, v0 big.Int
	//	v0.SetBit(&v0, from, 1)
	inf4 := [4]float64{math.Inf(1), math.Inf(1), math.Inf(1), math.Inf(1)}
	var score func(int, int, int)
	score = func(from, n, sx int) {
		if len(u.T[n]) == 1 {
			copy(cost4[n][:], inf4[:])
			cost4[n][u.Kmers[n][sx]>>1&3] = 0
			return
		}
		//fmt.Println("scoring node", n)
		for _, ch := range u.T[n] {
			if ch.To != from {
				score(n, ch.To, sx)
			}
		}
		for k := range "ACTG" {
			//fmt.Println("  base index", k)
			d := 0. // s(k)(v) from the text, the (dis)pars score for sym k
			for _, ch := range u.T[n] {
				if ch.To == from {
					continue
				}
				//fmt.Println("    child", ch, "dis[ch]", dis[ch])
				min := math.Inf(1)
				for i, dch := range cost4[ch.To] {
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

	// okay really label node n and compute cost for the edge leading to n.
	// axp is the "alphabet index of the parent"
	// return value is total cost (for the symbol position)
	var labelNodes func(int, graph.Half, int, int) float64
	labelNodes = func(from int, n graph.Half, sx, axp int) float64 {
		if len(u.T[n.To]) == 1 { // leaf nodes come with labels
			// just need to note distance
			if int(u.Kmers[n.To][sx]>>1&3) != axp {
				u.Costs[n.Label]++
				return 1
			}
			return 0
		}
		min := math.Inf(1)
		var xMin int
		for x, d := range cost4[n.To] {
			if x != axp {
				d++
			}
			if d < min {
				min = d
				xMin = x
			}
		}
		u.Kmers[n.To][sx] = "ACTG"[xMin]
		tot := 0.
		if axp >= 0 && xMin != axp {
			u.Costs[n.Label]++
			tot = 1
		}
		for _, ch := range u.T[n.To] {
			if ch.To != from {
				tot += labelNodes(n.To, ch, sx, xMin)
			}
		}
		return tot
	}

	// outer loop goes by symbol position
	tot := 0.
	for sx := range u.Kmers[0] {
		// depth first pass 1: accumulate dis-parsimony scores
		score(from, to.To, sx)
		score(to.To, from, sx)
		// compute (alphabet index of) most parsimonious symbol
		dMin := math.Inf(1)
		var axMin int
		for k := range "ACTG" {
			d := 0.
			for _, ch := range []int{from, to.To} {
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
			if d < dMin {
				dMin = d
				axMin = k
			}
		}
		// depth-first pass 2: build kmers for internal nodes, set edge costs
		tot += labelNodes(from, to, sx, axMin)
		tot += labelNodes(to.To, graph.Half{from, to.Label}, sx, axMin)
	}
	return tot
}

func (u *UTree) MaxParsimonyUnrooted() float64 {
	kmMin := make(Kmers, len(u.T))
	for n, k := range u.Kmers {
		if len(u.T[n]) > 1 {
			kmMin[n] = append(DNA8{}, k...)
		}
	}
	ctMin := make([]float64, len(u.Costs))
	minCost := math.Inf(1) // return value
	var df func(int, graph.Half)
	df = func(from int, n graph.Half) {
		c := u.MaxParsimonyEdge(from, n)
		if c < minCost {
			minCost = c
			for n, k := range u.Kmers {
				if len(u.T[n]) > 1 {
					copy(kmMin[n], k)
				}
			}
			copy(ctMin, u.Costs)
		}
		for _, to := range u.T[n.To] {
			if to.To != from {
				df(n.To, to)
			}
		}
	}
	for _, to := range u.T[0] {
		df(0, to)
	}
	for n, k := range kmMin {
		if len(u.T[n]) > 1 {
			copy(u.Kmers[n], k)
		}
	}
	copy(u.Costs, ctMin)
	return minCost
}

func (u *UTree) SwapEdges(a, ax, b, bx int) {
	aTo := u.T[a]
	bTo := u.T[b]

	// flip a's reciprocal arc to b
	to := u.T[aTo[ax].To] // neighbors of neighbor of a at index ax
	rx := 0
	for to[rx].To != a { // find index of reciprocal arc
		rx++
	}
	to[rx].To = b // flip

	// flip b's reciprocal arc to a
	to = u.T[bTo[bx].To]
	rx = 0
	for to[rx].To != b {
		rx++
	}
	to[rx].To = a

	// swap arcs from a and b
	aTo[ax], bTo[bx] = bTo[bx], aTo[ax]
}

func (u *UTree) SubTrees(a, b int) (w, x, y, z int) {
	// set w, x to neighbor indexes of a
	switch to := u.T[a]; b {
	case to[0].To:
		w, x = 1, 2
	case to[1].To:
		w, x = 0, 2
	default:
		w, x = 0, 1
	}
	// set y, z to neighbor indexes of b
	switch to := u.T[b]; a {
	case to[0].To:
		y, z = 1, 2
	case to[1].To:
		y, z = 0, 2
	default:
		y, z = 0, 1
	}
	return
}

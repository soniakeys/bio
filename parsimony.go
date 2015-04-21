package bio

import (
	"math/big"

	"github.com/soniakeys/graph"
)

var actg = DNA8("ACTG")

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
func DNA8MaxParsimonyRooted(tree graph.AdjacencyList, leaves Kmers) (kmers Kmers, costs []int) {
	root := len(tree) - 1
	m := len(leaves[0])
	kmers = append(leaves, make(Kmers, len(tree)-len(leaves))...)
	for i := len(leaves); i < len(tree); i++ {
		kmers[i] = make(DNA8, m)
	}
	costs = make([]int, len(tree))     // cost is edge length from parent
	cost4 := make([][4]int, len(tree)) // (dis)parsimony scores
	// TODO this const is ugly (max int doesn't quite work) inf as a special
	// large value usually a bad idea.
	const inf int = 1e9
	infa := [4]int{inf, inf, inf, inf}
	var score func(int, int)
	score = func(n, sx int) {
		if len(tree[n]) == 0 {
			copy(cost4[n][:], infa[:])
			cost4[n][kmers[n][sx]>>1&3] = 0
			return
		}
		//fmt.Println("scoring node", n)
		for _, ch := range tree[n] {
			score(ch, sx)
		}
		for k := range "ACTG" {
			//fmt.Println("  base index", k)
			d := 0 // s(k)(v) from the text, the (dis)pars score for sym k
			for _, ch := range tree[n] {
				//fmt.Println("    child", ch, "dis[ch]", dis[ch])
				min := inf
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
		min := inf
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
		kmers[n][sx] = actg[xMin]
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

// DNA8MaxParsimonyUnrooted solves kmers and costs for maximum parsimony
// (minimum total cost) in a rooted phylogenic tree.
//
// The argument tree is a undirected graph, nodes 0:len(leaves) must be leaf
// nodes.  It must be a tree and so have len(tree)-1 edges, labeled with
// numbers 0:len(tree)-1.
//
// Returned kmers will be the leaves argument with appended kmers labeling
// internal nodes and the root.  Returned cost slice is indexed by the
// graph edge labels.
func DNA8MaxParsimonyUnrooted(tree graph.LabeledAdjacencyList, leaves Kmers) (kmers Kmers, costs []int) {
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

	var minC int = 1e9
	var minTree = make([][]int, len(rt))
	var minKmers Kmers
	var minCosts []int
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
			c := 0
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
	costs = make([]int, len(tree)-1)
	for n, l := range minEdgeLabels {
		costs[l] += minCosts[n]
	}
	return minKmers[:len(tree)], costs
}

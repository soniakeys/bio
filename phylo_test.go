package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
	"github.com/soniakeys/graph"
)

func ExampleParseNewick() {
	nk := "(mouse)dog:7.1;"
	rt, _ := bio.ParseNewick(nk)
	fmt.Println("Newick:", nk)
	fmt.Println("Node  Name   HasWeight  Weight  Parent  Path-len-to-root  IsLeaf")
	p := rt.List.Paths
	for n, nn := range rt.Nodes {
		fmt.Printf("\u00a0%2d   %-5s  %9t %7.1f %5d %11d %12d\n",
			n, nn.Name, nn.HasWeight, nn.Weight,
			p[n].From, p[n].Len, rt.List.Leaves.Bit(graph.NI(n)))
	}
	fmt.Println("Max path len to root:", rt.List.MaxLen)
	// Output:
	// Newick: (mouse)dog:7.1;
	// Node  Name   HasWeight  Weight  Parent  Path-len-to-root  IsLeaf
	//   0   dog         true     7.1    -1           1            0
	//   1   mouse      false     0.0     0           2            1
	// Max path len to root: 2
}

func ExamplePhyloRootedTree_Newick() {
	pl, _ := bio.ParseNewick("(dog,((elephant:3,mouse:1.2),robot),cat);")
	rt := pl.RootedTree()

	fmt.Println(rt.Newick())
	// Output:
	// (dog,((elephant:3,mouse:1.2),robot),cat);
}

func ExamplePhyloList_NodeMap() {
	pl, _ := bio.ParseNewick("(dog,cat)fish;")
	m := pl.NodeMap()
	fmt.Println(m["cat"])
	fmt.Println(m["fish"])
	_, ok := m["robot"]
	fmt.Println(ok)
	// Output:
	// 2
	// 0
	// false
}

func ExamplePhyloList_PathLen() {
	pl, _ := bio.ParseNewick("(dog,((elephant,mouse),cat),robot);")
	//           .
	//          /|\
	//         / | \
	//      dog  .  robot
	//          / \
	//         .   cat
	//        / \
	// elephant  mouse
	m := pl.NodeMap()
	fmt.Println(pl.PathLen(m["cat"], m["mouse"]))
	// Output:
	// 3
}

func ExamplePhyloRootedTree_CharacterTable() {
	pl, _ := bio.ParseNewick("(dog,((elephant,mouse),cat),robot);")
	rt := pl.RootedTree()
	//              0
	//             /|\
	//            / | \
	//       dog 1  2  7 robot
	//             / \
	//            3   6 cat
	//           / \
	// elephant 4   5 mouse
	fmt.Println("76543210")
	fmt.Println("--------")
	for _, c := range rt.CharacterTable() {
		fmt.Printf("%0*b\n", len(rt.Tree.AdjacencyList), &c)
	}
	// Output:
	// 76543210
	// --------
	// 00111000
	// 01111100
}

func ExamplePhyloRootedTree_CharacterTable_leafRoot() {
	pl, _ := bio.ParseNewick("((dog,((elephant,mouse),cat)))robot;")
	rt := pl.RootedTree()
	//    robot 0
	//           \
	//            1
	//           / \
	//      dog 2   3
	//             / \
	//            4   7 cat
	//           / \
	// elephant 5   6 mouse
	fmt.Println("76543210")
	fmt.Println("--------")
	for _, c := range rt.CharacterTable() {
		fmt.Printf("%0*b\n", len(rt.Tree.AdjacencyList), &c)
	}
	// Output:
	// 76543210
	// --------
	// 01110000
	// 11111000
}

func ExampleStrKmers_CharacterTable() {
	sk := bio.StrKmers{
		"GCATTACC",
		"TTCGTACC",
		"TCATGACC",
		"TCAGTCCC",
		"TCCGTATC",
	}
	fmt.Println("43210  Position")
	fmt.Println("-----  -")
	ct, pos, _ := sk.CharacterTable()
	for i, c := range ct {
		fmt.Printf("%0*b  %d\n", len(sk), &c, pos[i])
	}
	// Output:
	// 43210  Position
	// -----  -
	// 10010  2
	// 00101  3
}

func ExamplePhyloRootedTree_MaxParsimony() {
	//       6
	//     /   \
	//   4       5
	//  / \     / \
	// 0   1   2   3
	// |   |   |   ATGGACGA
	// |   |   CTGCGCTG
	// |   ATTGCGAC
	// CAAATCCC
	tree := graph.Directed{graph.AdjacencyList{
		4: {0, 1},
		5: {2, 3},
		6: {4, 5},
	}}
	nodes := make([]bio.PhyloRootedNode, len(tree.AdjacencyList))
	rt := &bio.PhyloRootedTree{
		Tree:  tree,
		Nodes: nodes,
		Root:  6,
	}
	kmers := bio.Kmers{
		bio.DNA8("CAAATCCC"),
		bio.DNA8("ATTGCGAC"),
		bio.DNA8("CTGCGCTG"),
		bio.DNA8("ATGGACGA"),
		bio.DNA8("--------"),
		bio.DNA8("--------"),
		bio.DNA8("--------"),
	}
	rt.MaxParsimony(kmers)
	treeWt := 0.
	fmt.Println("node   kmer    weight from parent")
	for i, k := range kmers {
		fmt.Printf("%d    %s         %g\n", i, k, nodes[i].Weight)
		treeWt += nodes[i].Weight
	}
	fmt.Println(">             total:", treeWt)
	// Output:
	// node   kmer    weight from parent
	// 0    CAAATCCC         5
	// 1    ATTGCGAC         3
	// 2    CTGCGCTG         5
	// 3    ATGGACGA         1
	// 4    ATAGACAC         1
	// 5    ATGGACAA         1
	// 6    ATAGACAA         0
	// >             total: 16
}

func ExamplePhyloList_DistanceMatrix() {
	l, _ := bio.ParseNewick("((:11,:2):4,:6,:7);")
	// gives:
	//         ----0----
	//  (wt=4)/    |    \
	//       1     |(6)  |(7)
	//  (11)/ \(2) |     |
	//     2   3   4     5
	m := l.DistanceMatrix([]graph.NI{2, 3, 4, 5})
	for _, r := range m {
		fmt.Printf("%2.0f\n", r)
	}
	fmt.Println("-------")
	m = l.DistanceMatrix([]graph.NI{5, 2})
	for _, r := range m {
		fmt.Printf("%2.0f\n", r)
	}
	// Output:
	// [ 0 13 21 22]
	// [13  0 12 13]
	// [21 12  0 13]
	// [22 13 13  0]
	// -------
	// [ 0 22]
	// [22  0]
}

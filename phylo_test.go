package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleParseNewick() {
	nk := "(mouse)dog:7.1;"
	rt, _ := bio.ParseNewick(nk)
	fmt.Println("Newick:", nk)
	fmt.Println("Node  Name   HasWeight  Weight  Parent  Path-len-to-root  IsLeaf")
	p := rt.Tree.Paths
	for n, nn := range rt.Nodes {
		fmt.Printf("\u00a0%2d   %-5s  %9t %7.1f %5d %11d %12d\n",
			n, nn.Name, nn.HasWeight, nn.Weight,
			p[n].From, p[n].Len, rt.Tree.Leaves.Bit(n))
	}
	fmt.Println("Max path len to root:", rt.Tree.MaxLen)
	// Output:
	// Newick: (mouse)dog:7.1;
	// Node  Name   HasWeight  Weight  Parent  Path-len-to-root  IsLeaf
	//   0   dog         true     7.1    -1           1            0
	//   1   mouse      false     0.0     0           2            1
	// Max path len to root: 2
}

func ExamplePhyloRootedTree_Newick() {
	pt, _ := bio.ParseNewick("(dog,((elephant:3,mouse:1.2),robot),cat);")
	rt := pt.RootedTree()

	fmt.Println(rt.Newick())
	// Output:
	// (dog,((elephant:3,mouse:1.2),robot),cat);
}

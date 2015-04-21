package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
	"github.com/soniakeys/graph"
)

func ExampleDNA8MaxParsimonyRooted() {
	tree := graph.AdjacencyList{
		4: []int{0, 1},
		5: []int{2, 3},
		6: []int{4, 5},
	}
	leaves := bio.Kmers{
		bio.DNA8("CAAATCCC"),
		bio.DNA8("ATTGCGAC"),
		bio.DNA8("CTGCGCTG"),
		bio.DNA8("ATGGACGA"),
	}
	kmers, costs := bio.DNA8MaxParsimonyRooted(tree, leaves)
	c := 0
	fmt.Println("node   label   cost from parent")
	for i, k := range kmers {
		fmt.Printf("%d    %s         %d\n", i, k, costs[i])
		c += costs[i]
	}
	fmt.Println(">             total:", c)
	// Output:
	// node   label   cost from parent
	// 0    CAAATCCC         5
	// 1    ATTGCGAC         3
	// 2    CTGCGCTG         5
	// 3    ATGGACGA         1
	// 4    ATAGACAC         1
	// 5    ATGGACAA         1
	// 6    ATAGACAA         0
	// >             total: 16
}

func ExampleDNA8MaxParsimonyUnrooted() {
	// diagram:
	//
	//   TCGGCCAA (node 0)                 (node 1) CCTGGCTG
	//              \                           /
	//               \ edge 0                  / edge 1
	//                \           edge 4      /
	//               (node 4)--------------(node 5)
	//                /                       \
	//               / edge 2                  \ edge 3
	//              /                           \
	//   CACAGGAT (node 2)                 (node 3) TGAGTACC
	//
	tree := graph.LabeledAdjacencyList{
		0: {{To: 4, Label: 0}},
		1: {{To: 4, Label: 1}},
		2: {{To: 5, Label: 2}},
		3: {{To: 5, Label: 3}},
		4: {{To: 0, Label: 0}, {To: 1, Label: 1}, {To: 5, Label: 4}},
		5: {{To: 2, Label: 2}, {To: 3, Label: 3}, {To: 4, Label: 4}},
	}
	leaves := bio.Kmers{
		bio.DNA8("TCGGCCAA"),
		bio.DNA8("CCTGGCTG"),
		bio.DNA8("CACAGGAT"),
		bio.DNA8("TGAGTACC"),
	}
	kmers, costs := bio.DNA8MaxParsimonyUnrooted(tree, leaves)
	fmt.Println("node  label")
	for n, k := range kmers {
		fmt.Println(n, "   ", k)
	}
	fmt.Println("edge  cost")
	for e, c := range costs {
		fmt.Println(e, "   ", c)
	}
	// Output:
	// node  label
	// 0     TCGGCCAA
	// 1     CCTGGCTG
	// 2     CACAGGAT
	// 3     TGAGTACC
	// 4     CCAGGCAA
	// 5     CAAGGAAA
	// edge  cost
	// 0     3
	// 1     3
	// 2     4
	// 3     5
	// 4     2
}

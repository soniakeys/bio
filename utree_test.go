package bio_test

import (
	"fmt"
	"reflect"
	"testing"

	"github.com/soniakeys/bio"
	"github.com/soniakeys/graph"
)

func ExampleDNA8MaxParsimonyRooted() {
	tree := graph.Directed{graph.AdjacencyList{
		4: {0, 1},
		5: {2, 3},
		6: {4, 5},
	}}
	leaves := bio.Kmers{
		bio.DNA8("CAAATCCC"),
		bio.DNA8("ATTGCGAC"),
		bio.DNA8("CTGCGCTG"),
		bio.DNA8("ATGGACGA"),
	}
	kmers, costs := bio.DNA8MaxParsimonyRooted(tree, leaves)
	c := 0.
	fmt.Println("node   label   cost from parent")
	for i, k := range kmers {
		fmt.Printf("%d    %s         %g\n", i, k, costs[i])
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

func ExampleUTree_MaxParsimonyUnrooted() {
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
	u := bio.UTree{
		NLeaves: 4,
		T: graph.LabeledUndirected{graph.LabeledAdjacencyList{
			0: {{To: 4, Label: 0}},
			1: {{To: 4, Label: 1}},
			2: {{To: 5, Label: 2}},
			3: {{To: 5, Label: 3}},
			4: {{To: 0, Label: 0}, {To: 1, Label: 1}, {To: 5, Label: 4}},
			5: {{To: 2, Label: 2}, {To: 3, Label: 3}, {To: 4, Label: 4}},
		}},
		Kmers: bio.Kmers{
			bio.DNA8("TCGGCCAA"),
			bio.DNA8("CCTGGCTG"),
			bio.DNA8("CACAGGAT"),
			bio.DNA8("TGAGTACC"),
			bio.DNA8("........"),
			bio.DNA8("........"),
		},
		Costs: make([]float64, 5),
	}
	u.MaxParsimonyUnrooted()
	fmt.Println("node  label")
	for n, k := range u.Kmers {
		fmt.Println(n, "   ", k)
	}
	fmt.Println("edge  cost")
	for e, c := range u.Costs {
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

func TestSwapEdges(t *testing.T) {
	//     (10)                 (11)
	//      \                   /
	//       \ 20              / 21
	//        \       22      /
	//        (12)------------(13)
	//        /               \
	//       / 23              \ 24
	//      /                   \
	//     (14)                  (15)
	//      \                   /
	//       \ 25              / 26
	//        \               /
	//        (16)            (17)
	//
	// Test will swap lower edges from nodes 12 and 13.
	u := bio.UTree{T: graph.LabeledUndirected{graph.LabeledAdjacencyList{
		10: {{12, 20}},
		11: {{13, 21}},
		12: {{10, 20}, {13, 22}, {14, 23}},
		13: {{11, 21}, {15, 24}, {12, 22}},
		14: {{12, 23}, {16, 25}},
		15: {{13, 24}, {17, 26}},
		16: {{14, 25}},
		17: {{15, 26}},
	}}}
	u.SwapEdges(
		12, // node 12
		2,  // the lower edge is to node 14, at neighbor index 2
		13, // node 13
		1,  // lower edge is to node 15, at index 1
	)
	want := graph.LabeledUndirected{graph.LabeledAdjacencyList{
		10: {{12, 20}},
		11: {{13, 21}},
		// 12 should now have a half edge to 15, labeled 24
		12: {{10, 20}, {13, 22}, {15, 24}},
		// 13 should now have a half edge to 14, labeled 23
		13: {{11, 21}, {14, 23}, {12, 22}},
		// 14 should now have a half edge back to 13, labeled 23
		14: {{13, 23}, {16, 25}},
		// 15 should now have a half edge back to 12, labeled 24
		15: {{12, 24}, {17, 26}},
		16: {{14, 25}},
		17: {{15, 26}},
	}}
	if !reflect.DeepEqual(u.T, want) {
		t.Fatal("rats.")
	}
	// large parsimony doesn't really care how the edge labels are swapped
	// since it's recomputing costs anyway.  an alternative function would
	// be swapSubTrees, where the edge labels stayed with a and b and just
	// the nodes at the ends of the edges were swapped.  it's the same amount
	// of work either way though, and swapping "closer" to a and b seems more
	// intuitive.  either way, edge labels do need to be left in a consistent
	// state, so this tested behavior is good.
}

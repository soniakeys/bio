package bio_test

import (
	"fmt"
	"log"

	"github.com/soniakeys/bio"
	"github.com/soniakeys/graph"
)

var strSDB = bio.Str("TAATTATTAA")

func ExampleStr_DeBruijn() {
	fmt.Println("str:", strSDB)
	g, jmers := strSDB.DeBruijn(4)
	fmt.Println("DeBruijn graph:")
	for n, to := range g {
		fmt.Println(jmers[n], n, to)
	}
	// an Eulerian path exists in a DeBruijn graph
	p, err := graph.AdjacencyList(g).EulerianPath()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("Reconstructed path:", p)
	r, err := jmers.OverlapKmers(p)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("Reconstructed str:", r)
	// Output:
	// str: TAATTATTAA
	// DeBruijn graph:
	// TAA 0 [1]
	// AAT 1 [2]
	// ATT 2 [3 3]
	// TTA 3 [4 0]
	// TAT 4 [2]
	// Reconstructed path: [0 1 2 3 4 2 3 0]
	// Reconstructed str: TAATTATTAA
}

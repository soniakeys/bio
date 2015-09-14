package bio_test

import (
	"fmt"
	"log"

	"github.com/soniakeys/bio"
	"github.com/soniakeys/graph"
)

var strKDB = bio.Str("TAATTATTAA")

func ExampleStrKmers_DeBruijn() {
	fmt.Println("str:", strKDB)
	// extract kmers
	k := 4
	fmt.Printf("%d-mers:\n", k)
	kmers := make(bio.StrKmers, len(strKDB)-k+1)
	for i, j := 0, k; j <= len(strKDB); i, j = i+1, j+1 {
		kmers[i] = strKDB[i:j]
		fmt.Println(kmers[i])
	}

	g, jmers, err := kmers.DeBruijn()
	if err != nil {
		log.Fatal(err)
	}
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
	// 4-mers:
	// TAAT
	// AATT
	// ATTA
	// TTAT
	// TATT
	// ATTA
	// TTAA
	// DeBruijn graph:
	// TAA 0 [1]
	// AAT 1 [2]
	// ATT 2 [3 3]
	// TTA 3 [4 0]
	// TAT 4 [2]
	// Reconstructed path: [0 1 2 3 4 2 3 0]
	// Reconstructed str: TAATTATTAA
}

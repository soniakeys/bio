package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
	"github.com/soniakeys/multiset"
)

func ExampleSolvePartialDigest() {
	l := multiset.Multiset{}
	l.AddElements(2, 2, 3, 3, 4, 5, 6, 7, 8, 10)
	s := bio.SolvePartialDigest(l)
	for _, s1 := range s {
		fmt.Println(s1)
	}
	n := len(s[0])
	fmt.Println(n)
	fmt.Println(n*(n-1)/2, l.Cardinality())
	// Output:
	// [0 10 8 3 6]
	// [0 10 2 7 4]
	// 5
	// 10 10
}

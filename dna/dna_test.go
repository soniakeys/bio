package dna_test

import (
	"fmt"

	"github.com/soniakeys/bio/dna"
)

func ExampleSeq_Freq() {
	s := dna.Seq("Agatha")
	h := s.Freq()
	fmt.Println(len(h), "symbols")
	for _, b := range []byte("Aaght") {
		fmt.Printf("%c %d\n", b, h[b])
	}
	// Output:
	// 5 symbols
	// A 1
	// a 2
	// g 1
	// h 1
	// t 1
}

func ExampleStrict_BaseFreq() {
	s := dna.Strict("Gattaca")
	a, c, t, g := s.BaseFreq()
	fmt.Println("a", a)
	fmt.Println("c", c)
	fmt.Println("t", t)
	fmt.Println("g", g)
	// Output:
	// a 3
	// c 1
	// t 2
	// g 1
}

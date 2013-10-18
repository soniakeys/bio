package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleDNAStrict_BaseFreq() {
	s := bio.DNAStrict("Gattaca")
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

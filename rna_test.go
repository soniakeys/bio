package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleRNAStrict_BaseFreq() {
	s := bio.RNAStrict("Agua")
	a, c, u, g := s.BaseFreq()
	fmt.Println("a", a)
	fmt.Println("c", c)
	fmt.Println("u", u)
	fmt.Println("g", g)
	// Output:
	// a 2
	// c 0
	// u 1
	// g 1
}

func ExampleRNAStrict_Translate() {
	s := bio.RNAStrict("AugGcgAacAauUacUga")
	a, err := s.Translate()
	fmt.Printf("%T\n", a)
	fmt.Println(a)
	if err != nil {
		fmt.Println(err)
	}
	// Output:
	// bio.AAStrict
	// MANNY
}

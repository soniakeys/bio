package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleFreq() {
	s := bio.DNA("Agatha")
	h := bio.Freq(s)
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

func ExampleHamming() {
	s := bio.DNA("AGGCTTAC")
	t := bio.DNA("AGgCTAAC")
	h, err := bio.Hamming(s, t)
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Println(h)
	// Output:
	// 2
}

func ExampleAllIndex() {
	s := bio.DNA("Atatgatatat")
	m := bio.DNA("atat")
	fmt.Println(bio.AllIndex(s, m))
	// Output:
	// [6 8]
}

func ExampleToUpper() {
	t := bio.DNA("AGgCT-AC?")
	fmt.Println(bio.DNA(bio.ToUpper(t)))
	// Output:
	// AGGCT-AC?
}

func ExampleToLower() {
	t := bio.DNA("AGgCT-AC?")
	fmt.Println(bio.DNA(bio.ToLower(t)))
	// Output:
	// aggct-ac?
}

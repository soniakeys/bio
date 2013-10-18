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

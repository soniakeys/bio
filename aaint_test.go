package bio_test

import (
	"fmt"
	"sort"

	"github.com/soniakeys/bio"
)

func ExampleCut() {
	s := sort.IntSlice{3, 1, 4, 3, 2}
	fmt.Println(s)
	for i := 0; i <= 6; i++ {
		c := bio.Cut(s, i)
		fmt.Printf("first %d: %v", i, s[:c])
		switch {
		case c > i:
			fmt.Printf(" (%d way tie for last)\n", c-i+1)
		case c < i:
			fmt.Printf(" (only %d in list)\n", c)
		default:
			fmt.Println()
		}
	}
	// Output:
	// [3 1 4 3 2]
	// first 0: []
	// first 1: [1]
	// first 2: [1 2]
	// first 3: [1 2 3 3] (2 way tie for last)
	// first 4: [1 2 3 3]
	// first 5: [1 2 3 3 4]
	// first 6: [1 2 3 3 4] (only 5 in list)
}

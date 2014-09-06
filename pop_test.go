package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleLucasU() {
	// U(1,-1) are Fibonacci numbers
	for n := 0; n <= 6; n++ {
		fmt.Println(n, bio.LucasU(n, 1, -1))
	}
	// Output:
	// 0 0
	// 1 1
	// 2 1
	// 3 2
	// 4 3
	// 5 5
	// 6 8
}

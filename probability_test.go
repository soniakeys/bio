package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleMendel1() {
	d, h, r := bio.Mendel1(25, 22, 27)
	fmt.Printf("%.5f %.5f %.5f\n", d, h, r)
	// Output:
	// 0.23427 0.50444 0.26129
}

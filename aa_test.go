package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleAAInt_LinearSpec() {
	fmt.Println(bio.AA20("NQEL").Int().LinearSpec())
	// Output:
	// [0 113 114 128 129 242 242 257 370 371 484]
}

package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleBWT_AllIndex() {
	b := bio.NewBWT("AATCGGGTTCAATCGGGGT", 0, 5)
	fmt.Println(b.AllIndex("ATCG"))
	fmt.Println(b.AllIndex("GGGT"))
	// Output:
	// [11 1]
	// [15 4]
}

package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleSuffixArray() {
	sa := &bio.SuffixArray{}
	sa.AddSeq(bio.DNA8("GATTACA"))
	sa.AddSeq(bio.DNA8("TAGACCA"))
	sa.AddSeq(bio.DNA8("ATACA"))
	fmt.Println(string(sa.LongestCommonSubSeq()))
	// Output:
	// AC
}

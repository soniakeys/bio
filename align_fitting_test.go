package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

type matchAligner struct{}

func (matchAligner) Score(a, b byte) int {
	if a == b {
		return 1
	}
	return -1
}

func ExampleAlignPair_fitting() {
	s1 := bio.Seq("TAGATA")
	s2 := bio.Seq("GTAGGCTTAAGGTTA")

	score, t1, t2 := bio.AlignPair("fitting", s1, s2, matchAligner{}, 1)
	fmt.Println(score)
	fmt.Println(t1)
	fmt.Println(t2)
	// Output:
	// 2
	// TAGA-T-A
	// TAGGCTTA
}

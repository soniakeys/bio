package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

type matchAligner2 struct {
	match    int
	mismatch int
}

func (m matchAligner2) Score(a, b byte) int {
	if a == b {
		return m.match
	}
	return m.mismatch
}

func ExampleAlignPair_overlap() {
	s1 := bio.Seq("PAWHEAE")
	s2 := bio.Seq("HEAGAWGHEE")

	score, t1, t2 := bio.AlignPair("overlap", s1, s2, matchAligner2{1, -2}, 2)
	fmt.Println(score)
	fmt.Println(t1)
	fmt.Println(t2)
	// Output:
	// 1
	// HEAE
	// HEA-
}

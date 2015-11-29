package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleLinearGap() {
	fmt.Println(bio.LinearGap(
		bio.Seq("PLEASANTLY"),
		bio.Seq("-MEA--N-LY"), bio.Blosum62, 5))
	fmt.Println(bio.LinearGap(
		bio.Seq("LYPRTEINSTRIN"),
		bio.Seq("LY---EINSTEIN"), bio.Pam250, 5))
	// Output:
	// 8
	// 23
}

func ExampleConstantGap() {
	fmt.Println(bio.ConstantGap(
		bio.Seq("PLEASANTLY"),
		bio.Seq("-MEA--N-LY"), bio.Blosum62, 5))
	// Output:
	// 13
}

func ExampleAlignGlobal() {
	s1 := bio.Seq("PLEASANTLY")
	s2 := bio.Seq("MEANLY")
	score, t1, t2 := bio.AlignGlobal(s1, s2, bio.Blosum62, 5)
	fmt.Println(score)
	fmt.Println(t1)
	fmt.Println(t2)
	// Output:
	// 8
	// PLEASANTLY
	// -MEA--N-LY
}

func ExampleAlignPair_global() {
	s1 := bio.Seq("PLEASANTLY")
	s2 := bio.Seq("MEANLY")
	score, t1, t2 := bio.AlignPair("global", s1, s2, bio.Blosum62, 5)
	fmt.Println(score)
	fmt.Println(t1)
	fmt.Println(t2)
	// Output:
	// 8
	// PLEASANTLY
	// -MEA--N-LY
}

func FAILING_ExampleAlignGlobalLS() {
	s1 := bio.Seq("PLEASANTLY")
	s2 := bio.Seq("MEANLY")
	score, t1, t2 := bio.AlignGlobalLS(s1, s2, bio.Blosum62, 5)
	fmt.Println(score)
	fmt.Println(t1)
	fmt.Println(t2)
	// Output:
	// 8
	// PLEASANTLY
	// -MEA--N-LY
}

func ExampleAlignLocal() {
	s1 := bio.Seq("MEANLY")
	s2 := bio.Seq("PENALTY")
	score, t1, t2 := bio.AlignLocal(s1, s2, bio.Pam250, 5)
	fmt.Println(score)
	fmt.Println(t1)
	fmt.Println(t2)

	s1 = bio.Seq("MEANLYPRTEINSTRING")
	s2 = bio.Seq("PLEASANTLYEINSTEIN")
	score, t1, t2 = bio.AlignLocal(s1, s2, bio.Pam250, 5)
	fmt.Println(score)
	fmt.Println(t1)
	fmt.Println(t2)

	// Output:
	// 15
	// EANL-Y
	// ENALTY
	// 23
	// LYPRTEINSTRIN
	// LY---EINSTEIN
}

func ExampleAlignPair_local() {
	s1 := bio.Seq("MEANLY")
	s2 := bio.Seq("PENALTY")
	score, t1, t2 := bio.AlignPair("local", s1, s2, bio.Pam250, 5)
	fmt.Println(score)
	fmt.Println(t1)
	fmt.Println(t2)

	s1 = bio.Seq("MEANLYPRTEINSTRING")
	s2 = bio.Seq("PLEASANTLYEINSTEIN")
	score, t1, t2 = bio.AlignPair("local", s1, s2, bio.Pam250, 5)
	fmt.Println(score)
	fmt.Println(t1)
	fmt.Println(t2)

	// Output:
	// 15
	// EANL-Y
	// ENALTY
	// 23
	// LYPRTEINSTRIN
	// LY---EINSTEIN
}

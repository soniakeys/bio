package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleKMP_Index() {
	k := bio.NewKMP(bio.Seq("ababaa"))
	fmt.Println(k.Index(bio.Seq("abaababaab")))
	// Output:
	// 3
}

func ExampleKMP_AllIndex() {
	k := bio.NewKMP(bio.Seq("ababaa"))
	fmt.Println(k.AllIndex(bio.Seq("abaababaababaa")))
	// Output:
	// [3 8]
}

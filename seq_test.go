package bio_test

import (
	"fmt"
	"sort"

	"github.com/soniakeys/bio"
)

func ExampleSeq_FreqMap() {
	s := bio.Seq("Agatha")
	m := s.FreqMap()
	fmt.Println(len(m), "symbols")
	// some hoops because maps are unordered...
	var i []int
	for k := range m {
		i = append(i, int(k))
	}
	sort.Ints(i)
	for _, k := range i {
		fmt.Printf("%c %d\n", k, m[byte(k)])
	}
	// Output:
	// 5 symbols
	// A 1
	// a 2
	// g 1
	// h 1
	// t 1
}

func ExampleSeq_Hamming() {
	s := bio.Seq("AGGCTTAC")
	t := bio.Seq("AGgCTAAC")
	fmt.Println(s.Hamming(t))
	// Output:
	// 2
}

func ExampleSeq_PDistance() {
	s := bio.Seq("AGGCTTAC")
	t := bio.Seq("AGgCTAAC")
	fmt.Println(s.PDistance(t))
	// Output:
	// 0.25
}

func ExampleSeq_AllIndex() {
	s := bio.Seq("Atatgatatat")
	m := bio.Seq("atat")
	fmt.Println(s.AllIndex(m))
	// Output:
	// [5 7]
}

func ExampleSeq_AllCount() {
	s := bio.Seq("Atatgatatat")
	m := bio.Seq("atat")
	fmt.Println(s.AllCount(m))
	// Output:
	// 2
}

func ExampleSeq_ToUpper() {
	t := bio.Seq("AGgCT-AC?")
	fmt.Println(t.ToUpper())
	// Output:
	// AGGCT-AC?
}

func ExampleSeq_ToLower() {
	t := bio.Seq("AGgCT-AC?")
	fmt.Println(t.ToLower())
	// Output:
	// aggct-ac?
}

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

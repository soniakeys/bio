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

func ExampleKmerComposition_k2() {
	c := bio.KmerComposition(2, "ACGCGGCTCTGAAA")
	// (sort for uniform output)
	u := make([]string, len(c))
	i := 0
	for k := range c {
		u[i] = k
		i++
	}
	sort.Strings(u)
	for _, k := range u {
		fmt.Println(k, c[k])
	}
	// Output:
	// AA 2
	// AC 1
	// CG 2
	// CT 2
	// GA 1
	// GC 2
	// GG 1
	// TC 1
	// TG 1
}

func ExampleKmerComposition_k5() {
	c := bio.KmerComposition(5, "ATATATAG")
	// (sort for uniform output)
	u := make([]string, len(c))
	i := 0
	for k := range c {
		u[i] = k
		i++
	}
	sort.Strings(u)
	for _, k := range u {
		fmt.Println(k, c[k])
	}
	// Output:
	// ATATA 2
	// TATAG 1
	// TATAT 1
}

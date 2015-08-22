package bio_test

import (
	"fmt"
	"reflect"
	"sort"
	"testing"

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

func ExampleKCompositionDist() {
	fmt.Println(bio.KCompositionDist(3, "ATATATAG", "GATATA"))
	// Output:
	// 4
}

func ExampleKCompositionDistMat() {
	l := []string{
		"ATATATAG",
		"ATATATA",
		"GATATA",
	}
	d := bio.KCompositionDistMat(3, l)
	for i, s := range l {
		fmt.Printf("%-8s %.0f\n", s, d[i])
	}
	// Output:
	// ATATATAG [0 1 4]
	// ATATATA  [1 0 3]
	// GATATA   [4 3 0]
}

func TestKCompDist(t *testing.T) {
	l := []string{
		"CGGACACACAAAAAGAAAGAAGAATTTTTAGGATCTTTTGTGTGCGAATAACTATGAGGAAGATT",
		"CGGACACACAAAAAGAAAGAAGAATTTTTAGGATCTTTTGTGTGCGAATAACTATGAGGAAGATTAA",
		"CGGACACACAAAAAGAAAAAAGATTTTTTAAGACTTTTTGTGTACGAGTAACTATGAGGAAGATTAAC",
	}
	const k = 4
	d := bio.KCompositionDistMat(k, l)
	for i := 1; i < len(d); i++ {
		for j := range d[:i] {
			if d[i][j] != d[j][i] {
				t.Fatal("not symmetric")
			}
			if d[i][j] != float64(bio.KCompositionDist(k, l[i], l[j])) {
				t.Fatal("disagreement with KCompDist")
			}
		}
	}
	m := make(bio.DNA8List, len(l))
	for i, s := range l {
		m[i] = bio.DNA8(s)
	}
	e := m.KCompositionDistMat(k)
	if !reflect.DeepEqual(d, e) {
		t.Fatal("disagreement between KCompositionDistMat's")
	}
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

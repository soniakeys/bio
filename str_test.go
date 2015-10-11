package bio_test

import (
	"fmt"
	"log"
	"reflect"
	"sort"
	"testing"

	"github.com/soniakeys/bio"
	"github.com/soniakeys/graph"
)

func ExampleStr_DNA8HammingVariants() {
	for _, v := range bio.Str("Act").DNA8HammingVariants(1, nil) {
		fmt.Println(v)
	}
	// Output:
	// Act
	// Cct
	// Tct
	// Gct
	// Aat
	// Att
	// Agt
	// Aca
	// Acc
	// Acg
}

func ExampleStr_KmerComposition_k2() {
	c := bio.Str("ACGCGGCTCTGAAA").KmerComposition(2)
	// (sort for uniform output)
	u := make([]string, len(c))
	i := 0
	for k := range c {
		u[i] = string(k)
		i++
	}
	sort.Strings(u)
	for _, k := range u {
		fmt.Println(k, c[bio.Str(k)])
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

func ExampleStr_KmerComposition_k5() {
	c := bio.Str("ATATATAG").KmerComposition(5)
	// (sort for uniform output)
	u := make([]string, len(c))
	i := 0
	for k := range c {
		u[i] = string(k)
		i++
	}
	sort.Strings(u)
	for _, k := range u {
		fmt.Println(k, c[bio.Str(k)])
	}
	// Output:
	// ATATA 2
	// TATAG 1
	// TATAT 1
}

func ExampleStr_ReadPairComposition() {
	f := bio.Str("TAATGCCATGGGATGTT").ReadPairComposition(3, 2)
	// sort for predictable output:
	s := make([]string, len(f.Freq))
	i := 0
	for p, m := range f.Freq {
		s[i] = fmt.Sprint(p, m)
		i++
	}
	sort.Strings(s)
	fmt.Println("d:", f.D)
	for _, p := range s {
		fmt.Println(p)
	}
	// Output:
	// d: 2
	// {AAT CAT} 1
	// {ATG ATG} 2
	// {CAT GAT} 1
	// {CCA GGA} 1
	// {GCC GGG} 1
	// {GGG GTT} 1
	// {TAA CCA} 1
	// {TGC TGG} 1
	// {TGG TGT} 1
}

func ExampleStr_KCompositionDist() {
	fmt.Println(bio.Str("ATATATAG").KCompositionDist(3, bio.Str("GATATA")))
	// Output:
	// 4
}

func ExampleStrList_KCompositionDistMat() {
	l := bio.StrList{
		"ATATATAG",
		"ATATATA",
		"GATATA",
	}
	d := l.KCompositionDistMat(3)
	for i, s := range l {
		fmt.Printf("%-8s %.0f\n", s, d[i])
	}
	// Output:
	// ATATATAG [0 1 4]
	// ATATATA  [1 0 3]
	// GATATA   [4 3 0]
}

func TestStrList_KCompositionDistMat(t *testing.T) {
	l := bio.StrList{
		"CGGACACACAAAAAGAAAGAAGAATTTTTAGGATCTTTTGTGTGCGAATAACTATGAGGAAGATT",
		"CGGACACACAAAAAGAAAGAAGAATTTTTAGGATCTTTTGTGTGCGAATAACTATGAGGAAGATTAA",
		"CGGACACACAAAAAGAAAAAAGATTTTTTAAGACTTTTTGTGTACGAGTAACTATGAGGAAGATTAAC",
	}
	const k = 4
	d := l.KCompositionDistMat(k)
	for i := 1; i < len(d); i++ {
		for j := range d[:i] {
			if d[i][j] != d[j][i] {
				t.Fatal("not symmetric")
			}
			if d[i][j] != float64(l[i].KCompositionDist(k, l[j])) {
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

// compare to ExampleDeBruijns.  An example is hard here because of the
// composition map so it's a test here.
func TestStrFreq_DeBruijn(t *testing.T) {
	s := bio.Str("ACGCGGCTCTGAAA")
	k := 4
	c := s.KmerComposition(k)
	g, jmers, err := c.DeBruijn()
	if err != nil {
		t.Fatal(err)
	}
	p, err := graph.AdjacencyList(g).EulerianPath()
	if err != nil {
		t.Fatal(err)
	}
	r, err := jmers.OverlapKmers(p)
	if err != nil {
		t.Fatal(err)
	}
	if bio.Str(r) != s {
		t.Fatalf("want %s got %s", s, r)
	}
}

func ExampleStrKmers_Contigs() {
	kmers := bio.StrKmers{
		"ABCD",
		"BCDE",
		"CDEF",
		"DEFG",
		"EFGL",
		"FGLM",
		"GLMN",
		"LMNO",
		"BCDH",
		"CDHI",
		"DHIJ",
		"HIJK",
		"IJKL",
		"JKLM",
		"KLMN",
	}
	s, err := kmers.Contigs()
	if err != nil {
		log.Fatal(err)
	}
	for _, seq := range s {
		fmt.Println(seq)
	}
	// Output:
	// ABCD
	// BCDEFGLMN
	// BCDHIJKLMN
	// LMNO
}

func ExampleStrList_DistanceMatrix() {
	l := bio.StrList{
		"TTTCCATTTA",
		"GATTCATTTC",
		"TTTCCATTTT",
		"GTTCCATTTA",
	}
	f := func(a, b bio.Str) float64 { return float64(a.Hamming(b)) }
	m := l.DistanceMatrix(f)
	for _, r := range m {
		fmt.Printf("%2.0f\n", r)
	}
	// Output:
	// [ 0  4  1  1]
	// [ 4  0  4  3]
	// [ 1  4  0  2]
	// [ 1  3  2  0]
}

func ExampleStrKmers_ReadBreak() {
	l := bio.StrKmers{
		"ABCDEF",
		"CDEFGH",
	}
	m := l.ReadBreak(3)
	// sort for predictable output
	s := make([]string, len(m))
	i := 0
	for kmer := range m {
		s[i] = string(kmer)
		i++
	}
	sort.Strings(s)
	for _, kmer := range s {
		fmt.Println(kmer, m[bio.Str(kmer)])
	}
	// Output:
	// ABC 1
	// BCD 1
	// CDE 2
	// DEF 2
	// EFG 1
	// FGH 1
}

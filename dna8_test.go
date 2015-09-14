package bio_test

import (
	"bytes"
	"fmt"
	"reflect"
	"sort"
	"testing"

	"github.com/soniakeys/bio"
)

func ExampleDNA8_Cmp() {
	fmt.Println("Rule 1, sequence length:")
	for _, tc := range []struct{ s, t string }{
		{"a", "a"},
		{"a", "aa"},
		{"aa", "a"},
	} {
		fmt.Println(tc.s, tc.t, bio.DNA8(tc.s).Cmp(bio.DNA8(tc.t)))
	}
	fmt.Println("Rule 2, base order:")
	for _, tc := range []struct{ s, t string }{
		{"t", "g"}, // order is package-typical ACTG
		{"a", "C"}, // different bases compare case-insensitve.
		{"A", "c"},
	} {
		fmt.Println(tc.s, tc.t, bio.DNA8(tc.s).Cmp(bio.DNA8(tc.t)))
	}
	fmt.Println("Rule 3, case breaks ties:")
	for _, tc := range []struct{ s, t string }{
		{"a", "A"},
		{"A", "a"},
	} {
		fmt.Println(tc.s, tc.t, bio.DNA8(tc.s).Cmp(bio.DNA8(tc.t)))
	}
	// Output:
	// Rule 1, sequence length:
	// a a 0
	// a aa -1
	// aa a 1
	// Rule 2, base order:
	// t g -1
	// a C -1
	// A c -1
	// Rule 3, case breaks ties:
	// a A -1
	// A a 1
}

func BenchmarkHammingVariantsRef(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariantsRef(2)
	}
}

func BenchmarkHammingVariants(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariants(2)
	}
}

func ExampleDNA8_FreqArray_k2() {
	fmt.Println(bio.DNA8("ACGCGGCTCTGAAA").FreqArray(2))
	// Output:
	// [2 1 0 0 0 0 2 2 0 1 0 1 1 2 0 1]
}

func ExampleDNA8_FreqArray_k5() {
	f := bio.DNA8("ATATATAG").FreqArray(5)
	fmt.Println(len(f))
	// Output:
	// 1024
}

func ExampleDNA8_ModalKmers() {
	s := bio.DNA8("ACGTTGCATGTCGCATGATGCATGAGAGCT")
	fmt.Println(s.ModalKmers(4))
	// Output:
	// [GCAT CATG]
}

func ExampleDNA8_ModalSmallKmers() {
	s := bio.DNA8("ACGTTGCATGTCGCATGATGCATGAGAGCT")
	fmt.Println(s.ModalSmallKmers(4))
	// Output:
	// [CATG GCAT]
}

func BenchmarkDNA8Complement100(b *testing.B) {
	for i := 0; i <= b.N; i++ {
		bio.DNA8Complement(d1k100[i%1000])
	}
}

func ExampleDNA8_AllIndex() {
	s := bio.DNA8("Atatgatatat")
	m := bio.DNA8("atat")
	fmt.Println(s.AllIndex(m))
	// Output:
	// [0 5 7]
}

func ExampleDNA8_AllCount() {
	s := bio.DNA8("Atatgatatat")
	m := bio.DNA8("atat")
	fmt.Println(s.AllCount(m))
	// Output:
	// 3
}

func ExampleDNA8_BaseFreq() {
	s := bio.DNA8("Gattaca")
	a, c, t, g := s.BaseFreq()
	fmt.Println("a", a)
	fmt.Println("c", c)
	fmt.Println("t", t)
	fmt.Println("g", g)
	// Output:
	// a 3
	// c 1
	// t 2
	// g 1
}

func ExampleDNA8_Transcribe() {
	d := bio.DNA8("Tact")
	r := d.Transcribe()
	fmt.Printf("%T %v\n", d, d)
	fmt.Printf("%T %v\n", r, r)
	// Output:
	// bio.DNA8 Tact
	// bio.RNA8 Uacu
}

func ExampleDNA8_ReverseComplement() {
	s := bio.DNA8("Atacaga")
	fmt.Println(s.ReverseComplement())
	// Output:
	// tctgtaT
}

func ExampleDNA8_GCContent() {
	s := bio.DNA8("cagggt")
	fmt.Printf("%.3f\n", s.GCContent())
	// Output:
	// 0.667
}

func ExampleGCSkew() {
	fmt.Println(bio.GCSkew('G'))
	fmt.Println(bio.GCSkew('c'))
	fmt.Println(bio.GCSkew('t'))
	// Output:
	// 1
	// -1
	// 0
}

func ExampleDNA8_MinGCSkew() {
	s := bio.DNA8("accagtgct")
	m := s.MinGCSkew()
	fmt.Println(m)
	w := bytes.Repeat([]byte{' '}, len(s))
	for _, x := range m {
		w[x] = '^'
	}
	fmt.Printf("|%s|\n", s)
	fmt.Printf("|%s|\n", w)
	// Output:
	// [2 3]
	// |accagtgct|
	// |  ^^     |
}

func ExampleKmers_Consensus() {
	c := bio.Kmers{
		bio.DNA8("GATTcca"),
		bio.DNA8("AATTcgg"),
		bio.DNA8("GACTaca"),
		bio.DNA8("GATAaca"),
		bio.DNA8("GATTaca"),
	}
	fmt.Println(c.Consensus())
	// Output:
	// GATTACA 28
}

func ExampleDNA8_PalFindAllIndex() {
	s := bio.DNA8("CAATGCATG")
	for _, p := range s.PalFindAllIndex(4, 8) {
		fmt.Printf("%d: %s\n", p.Index, s[p.Index:p.Index+p.Len])
	}
	// Output:
	// 3: TGCA
	// 2: ATGCAT
	// 5: CATG
}

func ExampleDNA8_HammingVariants() {
	seq := bio.DNA8("Act")
	vars := seq.HammingVariants(1)
	for _, v := range vars {
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

func ExampleNumHammingVariants() {
	fmt.Println(bio.NumHammingVariants(3, 1))
	fmt.Println(bio.NumHammingVariants(9, 3))
	fmt.Println(bio.NumHammingVariants(1000, 10))
	// Output:
	// 10
	// 2620
	// 15606547544311956375154416076
}

func ExampleDNA8_ModalHammingKmers() {
	s := bio.DNA8("aacaagcataaacattaaagag")
	for _, k := range s.ModalHammingKmers(5, 1) {
		fmt.Println(k)
	}
	// Output:
	// aaaaa
}

func ExampleDNA8_ModalHammingKmersRC() {
	s := bio.DNA8("aacaagcataaacattaaagag")
	for _, k := range s.ModalHammingKmersRC(5, 1) {
		fmt.Println(k)
	}
	// Output:
	// cttaa
	// ttaaa
	// tttaa
	// ttaag
	// aaaaa
	// ttttt
}

func ExampleDNA8_AAFindAllIndex() {
	s := bio.DNA8("gtgaaactttttccttggtttaatcaatat")
	p := bio.AA20("NQ")
	x := s.AAFindAllIndex(p)
	fmt.Println(x)

	y := x[0]
	fmt.Println(s[y : y+len(p)*3])
	y = x[1]
	fmt.Println(s[y : y+len(p)*3])
	// Output:
	// [21 14]
	// aatcaa
	// ttggtt
}

func ExampleDNA8_Hamming() {
	s1 := bio.DNA8("CATTAG")
	s2 := bio.DNA8("catagg")
	fmt.Println(s1.Hamming(s2))
	// Output:
	// 2
}

func ExampleDNA8_MotifHamming() {
	m := bio.DNA8("caa")
	s := bio.DNA8("gtgaaactt")
	fmt.Println(s.MotifHamming(m))
	// Output:
	// 1
}

func ExampleDNA8List_MotifHamming() {
	m := bio.DNA8("caa")
	l := bio.DNA8List{
		bio.DNA8("gtgaaactt"),
		bio.DNA8("gtgagactt"),
	}
	fmt.Println(l.MotifHamming(m))
	// Output:
	// 3
}

func ExampleDNA8_KmersNearestMotif() {
	m := bio.DNA8("caa")
	s := bio.DNA8("gtgaaactt")
	for _, k := range m.KmersNearestMotif(s) {
		fmt.Println(k)
	}
	// Output:
	// gaa
	// aaa
}

func ExampleDNA8_Inc() {
	m := bio.DNA8("ggt")
	fmt.Println("Initial motif:", m)
	fmt.Printf("Inc returns %5t, new motif: %s\n", m.Inc(), m)
	fmt.Printf("Inc returns %5t, new motif: %s\n", m.Inc(), m)
	fmt.Printf("Inc returns %5t, new motif: %s\n", m.Inc(), m)
	// Output:
	// Initial motif: ggt
	// Inc returns  true, new motif: ggg
	// Inc returns false, new motif: aaa
	// Inc returns  true, new motif: aac
}

func ExampleTiTvRatio8() {
	s := bio.DNA8("gggcttt")
	t := bio.DNA8("Atacata")
	fmt.Printf("%.3f\n", bio.TiTvRatio8(s, t))
	// Output:
	// 0.667
}

func TestHammingVariants(t *testing.T) {
	// for a few test cases, tests that HammingVariants result matches
	// HammingVariantsRef and test that len is as expected.
	tcs := []struct {
		kmer string
		d    int
	}{
		{"actg", 0},
		{"Act", 1},
		{"cGG", 2},
		{"atacaga", 3},
	}
	for _, tc := range tcs {
		ref := bio.DNA8(tc.kmer).HammingVariantsRef(tc.d)
		sort.Sort(bio.DNA8List(ref))
		got := bio.DNA8(tc.kmer).HammingVariants(tc.d)
		sort.Sort(bio.DNA8List(got))
		if !reflect.DeepEqual(ref, got) {
			t.Fatalf("tc %v", tc)
		}
		// test len
		want := bio.NumHammingVariants(len(tc.kmer), tc.d).Int64()
		if len(ref) != int(want) {
			t.Fatalf("len tc %v = %d, want %d", tc, len(ref), want)
		}
	}
}

func ExampleDNA8_UniqueKmers() {
	s := bio.DNA8("Atatatatag")
	u := s.UniqueKmers(3)
	// (sort for predictable output)
	o := make([]string, len(u))
	i := 0
	for k := range u {
		o[i] = string(k)
		i++
	}
	sort.Strings(o)
	for _, k := range o {
		fmt.Println(k)
	}
	// Output:
	// ATA
	// TAG
	// TAT
}

func ExampleDNA8_UniqueHammingKmers() {
	s := bio.DNA8("Atatatatag")
	u := s.UniqueHammingKmers(3, 1)
	// (sort for predictable output)
	o := make([]string, len(u))
	i := 0
	for k := range u {
		o[i] = string(k)
		i++
	}
	sort.Strings(o)
	for _, k := range o {
		fmt.Println(k)
	}
	// Output:
	// AAA
	// AAG
	// AAT
	// ACA
	// AGA
	// ATA
	// ATC
	// ATG
	// ATT
	// CAG
	// CAT
	// CTA
	// GAG
	// GAT
	// GTA
	// TAA
	// TAC
	// TAG
	// TAT
	// TCG
	// TCT
	// TGG
	// TGT
	// TTA
	// TTG
	// TTT
}

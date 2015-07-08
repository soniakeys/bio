package bio_test

import (
	"bytes"
	"fmt"
	"testing"

	"github.com/soniakeys/bio"
)

func BenchmarkHammingVariants1(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariants(1)
	}
}

func BenchmarkVarNeighbors1(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariants1(1)
	}
}

func BenchmarkHammingVariants2(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariants(2)
	}
}

func BenchmarkVarNeighbors2(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariants1(2)
	}
}

func BenchmarkHammingVariants3(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariants(3)
	}
}

func BenchmarkVarNeighbors3(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariants1(3)
	}
}

func BenchmarkHammingVariants4(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariants(4)
	}
}

func BenchmarkVarNeighbors4(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariants1(4)
	}
}

func BenchmarkHammingVariants5(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariants(5)
	}
}

func BenchmarkVarNeighbors5(b *testing.B) {
	m := bio.DNA8("CATGTCGCA")
	for i := 0; i < b.N; i++ {
		m.HammingVariants1(5)
	}
}

func ExampleDNA8_FreqArray() {
	fmt.Println(bio.DNA8("ACGCGGCTCTGAAA").FreqArray(2))
	// Output:
	// [2 1 0 0 0 0 2 2 0 1 0 1 1 2 0 1]
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
	for _, v := range bio.DNA8("Act").HammingVariants(1) {
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

func ExampleDNA8_MotifSeqDist() {
	m := bio.DNA8("caa")
	s := bio.DNA8("gtgaaactt")
	fmt.Println(s.MotifHamming(m))
	// Output:
	// 1
}

func ExampleDNA8_MotifSetDist() {
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
	fmt.Println(m)
	m.Inc()
	fmt.Println(m)
	m.Inc()
	fmt.Println(m)
	m.Inc()
	fmt.Println(m)
	// Output:
	// ggt
	// ggg
	// aaa
	// aac
}

func ExampleTiTvRatio8() {
	s := bio.DNA8("gggcttt")
	t := bio.DNA8("Atacata")
	fmt.Printf("%.3f\n", bio.TiTvRatio8(s, t))
	// Output:
	// 0.667
}

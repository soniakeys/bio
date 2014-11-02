package bio_test

import (
	"bytes"
	"fmt"
	"math/rand"
	"testing"

	"github.com/soniakeys/bio"
)

var d1k90 []byte
var d1k100 []byte

func init() {
	d1k90 = make([]byte, 1000)
	d1k100 = make([]byte, 1000)
	// make 90% dna
	for i := range d1k90 {
		if rand.Intn(10) == 0 {
			d1k90[i] = byte(rand.Intn(256))
		} else {
			d1k90[i] = "actg"[rand.Intn(4)]
		}
	}
	// make 100% dna
	for i := range d1k100 {
		d1k100[i] = "actg"[rand.Intn(4)]
	}
}

func BenchmarkDNAComplement90(b *testing.B) {
	for i := 0; i <= b.N; i++ {
		bio.DNAComplement(d1k90[i%1000])
	}
}

func BenchmarkDNAComplement100(b *testing.B) {
	for i := 0; i <= b.N; i++ {
		bio.DNAComplement(d1k100[i%1000])
	}
}

func BenchmarkDNA8Complement100(b *testing.B) {
	for i := 0; i <= b.N; i++ {
		bio.DNA8Complement(d1k100[i%1000])
	}
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

func ExampleDNA_BaseFreq() {
	s := bio.DNA("Arithmetic")
	a, c, t, g := s.BaseFreq()
	fmt.Println("a", a)
	fmt.Println("c", c)
	fmt.Println("t", t)
	fmt.Println("g", g)
	// Output:
	// a 1
	// c 1
	// t 2
	// g 0
}

func ExampleDNA_Transcribe() {
	d := bio.DNA("Tristan")
	r := d.Transcribe()
	fmt.Printf("%T %v\n", d, d)
	fmt.Printf("%T %v\n", r, r)
	// Output:
	// bio.DNA Tristan
	// bio.RNA Urisuan
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

func ExampleDNA_ReverseComplement() {
	s := bio.DNA("Atacama")
	fmt.Println(s.ReverseComplement())
	// Output:
	// tktgtaT
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

func ExampleDNA_GCContent() {
	s := bio.DNA("carrot")
	fmt.Printf("%.3f\n", s.GCContent())
	// Output:
	// 0.333
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

func ExampleCountProfile() {
	set := []string{
		"GAT..ca",
		"AA##cgg",
		"GACrcca",
		"Gxxxaca",
		"GATmaca",
	}
	pm := make(bio.CountProfile, len(set[0]))
	for _, s := range set {
		pm.Add(bio.DNA(s))
	}
	fmt.Println(pm.Consensus())
	// Output:
	// GAT-ACA 20
}

func ExampleDNAConsensus() {
	c := []bio.DNA{
		bio.DNA("GAT..ca"),
		bio.DNA("AA##cgg"),
		bio.DNA("GACrcca"),
		bio.DNA("Gxxxaca"),
		bio.DNA("GATmaca"),
	}
	fmt.Println(bio.DNAConsensus(c))
	// Output:
	// GAT-ACA 20
}

func ExampleDNA8List_Consensus() {
	c := bio.DNA8List{
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
	m := s[y : y+len(p)*3]
	fmt.Println(m)
	t, _ := m.Translate()
	fmt.Println(t)
	// Output:
	// [21]
	// aatcaa
	// NQ
}

func ExampleDNA8_AAFindAllIndexRC() {
	s := bio.DNA8("gtgaaactttttccttggtttaatcaatat")
	p := bio.AA20("NQ")
	x := s.AAFindAllIndexRC(p)
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
	fmt.Println(m.MotifSeqDist(s))
	// Output:
	// 1
}

func ExampleDNA8_MotifSetDist() {
	m := bio.DNA8("caa")
	l := []bio.DNA8{
		bio.DNA8("gtgaaactt"),
		bio.DNA8("gtgagactt"),
	}
	fmt.Println(m.MotifSetDist(l))
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

func ExampleTiTvRatio() {
	s := bio.DNA("gggcttt")
	t := bio.DNA("Atacama")
	fmt.Printf("%.3f\n", bio.TiTvRatio(s, t))
	// Output:
	// 0.667
}

func ExampleTiTvRatio8() {
	s := bio.DNA8("gggcttt")
	t := bio.DNA8("Atacata")
	fmt.Printf("%.3f\n", bio.TiTvRatio8(s, t))
	// Output:
	// 0.667
}

package bio_test

import (
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

func ExampleDNA_ReverseComplement() {
	s := bio.DNA("Atacama")
	fmt.Println(s.ReverseComplement())
	// Output:
	// tktgtaT
}

func ExampleDNA_GCContent() {
	s := bio.DNA("carrot")
	fmt.Printf("%.3f\n", s.GCContent())
	// Output:
	// 0.333
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

func ExampleTiTvRatio() {
	s := bio.DNA("gggcttt")
	t := bio.DNA("Atacama")
	fmt.Printf("%.3f\n", bio.TiTvRatio(s, t))
	// Output:
	// 0.667
}

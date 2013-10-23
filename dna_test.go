package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleDNAStrict_BaseFreq() {
	s := bio.DNAStrict("Gattaca")
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

func ExampleDNA_Transcribe() {
	d := bio.DNA("Tristan")
	r := d.Transcribe()
	fmt.Printf("%T %v\n", d, d)
	fmt.Printf("%T %v\n", r, r)
	// Output:
	// bio.DNA Tristan
	// bio.RNA Urisuan
}

func ExampleDNAStrict_Transcribe() {
	d := bio.DNAStrict("Tact")
	r := d.Transcribe()
	fmt.Printf("%T %v\n", d, d)
	fmt.Printf("%T %v\n", r, r)
	// Output:
	// bio.DNAStrict Tact
	// bio.RNAStrict Uacu
}

func ExampleDNA_ReverseComplement() {
	s := bio.DNA("Atacama")
	fmt.Println(s.ReverseComplement())
	// Output:
	// tmtgtaT
}

func ExampleDNAStrict_ReverseComplement() {
	s := bio.DNAStrict("Atacaga")
	fmt.Println(s.ReverseComplement())
	// Output:
	// tctgtaT
}

func ExampleDNAStrict_CGFraction() {
	s := bio.DNAStrict("cagggt")
	fmt.Printf("%.3f\n", s.CGFraction())
	// Output:
	// 0.667
}

func ExampleDNAConsensus() {
	c := []bio.DNA{
		bio.DNA("GAT..CA"),
		bio.DNA("AA##CGG"),
		bio.DNA("GACrCCA"),
		bio.DNA("GxxxACA"),
		bio.DNA("GATmACA")}
	fmt.Println(bio.DNAConsensus(c))
	// Output:
	// GAT-ACA 20
}

func ExampleDNAStrictConsensus() {
	c := []bio.DNAStrict{
		bio.DNAStrict("GATTCCA"),
		bio.DNAStrict("AATTCGG"),
		bio.DNAStrict("GACTACA"),
		bio.DNAStrict("GATAACA"),
		bio.DNAStrict("GATTACA")}
	fmt.Println(bio.DNAStrictConsensus(c))
	// Output:
	// GATTACA 28
}

package bio_test

import (
	"fmt"

	"github.com/soniakeys/bio"
)

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
	// tmtgtaT
}

func ExampleDNA8_ReverseComplement() {
	s := bio.DNA8("Atacaga")
	fmt.Println(s.ReverseComplement())
	// Output:
	// tctgtaT
}

func ExampleDNA8_CGFraction() {
	s := bio.DNA8("cagggt")
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

func ExampleDNA8Consensus() {
	c := []bio.DNA8{
		bio.DNA8("GATTCCA"),
		bio.DNA8("AATTCGG"),
		bio.DNA8("GACTACA"),
		bio.DNA8("GATAACA"),
		bio.DNA8("GATTACA")}
	fmt.Println(bio.DNAConsensus8(c))
	// Output:
	// GATTACA 28
}

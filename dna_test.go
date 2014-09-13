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

func ExampleDNA8_GCContent() {
	s := bio.DNA8("cagggt")
	fmt.Printf("%.3f\n", s.GCContent())
	// Output:
	// 0.667
}

func ExampleDNAProfileMatrix() {
	set := []string{
		"GAT..ca",
		"AA##cgg",
		"GACrcca",
		"Gxxxaca",
		"GATmaca",
	}
	pm := make(bio.DNAProfileMatrix, len(set[0]))
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

func ExampleDNA8Consensus() {
	c := []bio.DNA8{
		bio.DNA8("GATTcca"),
		bio.DNA8("AATTcgg"),
		bio.DNA8("GACTaca"),
		bio.DNA8("GATAaca"),
		bio.DNA8("GATTaca"),
	}
	fmt.Println(bio.DNA8Consensus(c))
	// Output:
	// GATTACA 28
}

func ExampleDNA8_FindAllPalIndex() {
	s := bio.DNA8("CAATGCATG")
	for _, p := range s.FindAllPalIndex(4, 8) {
		fmt.Printf("%d: %s\n", p.Index, s[p.Index:p.Index+p.Len])
	}
	// Output:
	// 3: TGCA
	// 2: ATGCAT
	// 5: CATG
}

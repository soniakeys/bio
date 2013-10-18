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

package bio_test

import (
	"fmt"
	"sort"

	"github.com/soniakeys/bio"
)

func ExampleIsDNAStart() {
	fmt.Println(bio.IsDNAStart('A', 'G', 'G'))
	fmt.Println(bio.IsDNAStart('A', 'T', 'G'))
	fmt.Println(bio.IsDNAStart('a', 't', 'g'))
	// Output:
	// false
	// true
	// true
}

func ExampleIsRNAStart() {
	fmt.Println(bio.IsRNAStart('A', 'T', 'G'))
	fmt.Println(bio.IsRNAStart('A', 'U', 'G'))
	fmt.Println(bio.IsRNAStart('a', 'u', 'g'))
	// Output:
	// false
	// true
	// true
}

func ExampleRNA8_Translate() {
	s := bio.RNA8("AugGcgAacAauUacUga")
	a, err := s.Translate()
	fmt.Printf("%T\n", a)
	fmt.Println(a)
	if err != nil {
		fmt.Println(err)
	}
	// Output:
	// bio.AA20
	// MANNY
}

func ExampleDNA8_TranslateORF() {
	s := bio.DNA8("ctATGatcttctCATactacacATGtaCATaacccc")
	a := s.TranslateORF()
	// Sort result strings just so output order is repeatable.
	t := make([]string, len(a))
	for i, seq := range a {
		t[i] = seq.String()
	}
	sort.Strings(t)
	for _, s := range t {
		fmt.Println(s)
	}
	// Output:
	// MCSMRRS
	// MIFSYYTCT
	// MRRS
	// MYMCSMRRS
}

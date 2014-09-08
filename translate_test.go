package bio_test

import (
	"fmt"
	"sort"

	"github.com/soniakeys/bio"
)

func ExampleIsStartCodon() {
	fmt.Println(bio.IsStartCodon('C', 'U', 'G'))
	fmt.Println(bio.IsStartCodon('A', 'U', 'G'))
	fmt.Println(bio.IsStartCodon('a', 'u', 'g'))
	fmt.Println(bio.IsStartCodon('a', 't', 'g'))
	// Output:
	// false
	// true
	// true
	// true
}

func ExampleRNA8_Translate() {
	fmt.Println(bio.RNA8("AugGcgAacAauUacUga").Translate())
	// Output:
	// MANNY true
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

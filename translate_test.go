package bio_test

import (
	"fmt"
	"sort"

	"github.com/soniakeys/bio"
)

func ExampleRNAStrict_Translate() {
	s := bio.RNAStrict("AugGcgAacAauUacUga")
	a, err := s.Translate()
	fmt.Printf("%T\n", a)
	fmt.Println(a)
	if err != nil {
		fmt.Println(err)
	}
	// Output:
	// bio.AAStrict
	// MANNY
}

func ExampleDNAStrict_TranslateORF() {
	s := bio.DNAStrict("ctATGatcttctCATactacacATGtaCATaacccc")
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

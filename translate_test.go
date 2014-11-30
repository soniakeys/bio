package bio_test

import (
	"fmt"
	"io/ioutil"
	"sort"
	"testing"

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

func BenchmarkTranslateTable(b *testing.B) {
	d, err := ioutil.ReadFile("Vibrio_cholerae.txt")
	if err != nil {
		b.Skip("Vibrio_cholerae.txt not present")
	}
	seq := bio.DNA8(d)
	sn := len(seq) - 3 // (last byte in file is a \n)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p := i % sn
		_ = bio.TranslateCodon(seq[p], seq[p+1], seq[p+2])
	}
}

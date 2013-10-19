// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"errors"
)

// rna.go
//
// Types and methods that are specific to RNA or optimized for RNA.

// RNA type represents a RNA sequence.
//
// It is expected to hold RNA base symbols but other symbols are allowed.
// Methods on the type accomodate non-base symbols.
type RNA []byte

// RNAStrict type represents a sequence consisting strictly of RNA symbols.
//
// Allowed symbols are ACUGacug.  Methods on the type assume this.  Methods
// are case-insensitive but may produce nonsense results if the sequence
// contains non-base symbols.
type RNAStrict []byte

// String satisfies fmt.Stringer.
func (s RNA) String() string {
	return string(s)
}

// String satisfies fmt.Stringer.
func (s RNAStrict) String() string {
	return string(s)
}

// BaseFreq returns counts of each of the four RNA bases.
func (s RNAStrict) BaseFreq() (a, c, u, g int) {
	return baseFreq(s)
}

// lookup table for translation
const (
	aaStop = '*'
	codons = "KNNKTTTTIIIMRSSRQHHQPPPPLLLLRRRR*YY*SSSSLFFL*CCWEDDEAAAAVVVVGGGG"
)

func codonIndex(b0, b1, b2 byte) int {
	return int(b0&6<<3 | b1&6<<1 | b2&6>>1)
}

// RNAStart is the RNA start codon, AUG
var RNAStart = RNAStrict("AUG")

func IsRNAStart(s []byte) bool {
	return s[0]&0xdf == RNAStart[0] &&
		s[1]&0xdf == RNAStart[1] &&
		s[2]&0xdf == RNAStart[2]
}

// Translate translates RNA and returns an amino acid sequence.
//
// Translation begins at the start of the sequence.  The sequence should begin
// with the RNA start codon.  Translation ends at a stop codon or at the
// end of the sequence.
//
// Errors are returned for no start codon and no stop codon but the
// translated sequence is returned in any case.
func (s RNAStrict) Translate() (a AAStrict, err error) {
	if !IsRNAStart(s) {
		err = errors.New("No start codon")
	}
	for p := 0; ; p += 3 {
		if p+3 > len(s) {
			err = errors.New("No stop codon")
			return
		}
		aa := codons[codonIndex(s[p], s[p+1], s[p+2])]
		if aa == aaStop {
			return
		}
		a = append(a, aa)
	}
}

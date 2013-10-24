// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

// translate.go
//
// Definitions concerning tranlation of DNA and RNA into amino acid sequences.

package bio

import (
	"errors"
	"regexp"
)

// lookup table for translation.
// condonIndex and translateCodon should inline nicely.
const (
	aaStop = '*'
	codons = "KNNKTTTTIIIMRSSRQHHQPPPPLLLLRRRR*YY*SSSSLFFL*CCWEDDEAAAAVVVVGGGG"
)

func codonIndex8(b0, b1, b2 byte) int {
	return int(b0&6<<3 | b1&6<<1 | b2&6>>1)
}

func translateCodon8(b0, b1, b2 byte) byte {
	return codons[codonIndex8(b0, b1, b2)]
}

// RNAStart is the RNA start codon, AUG
const RNAStart = "AUG"

// DNAStart is the DNA start codon, ATG
const DNAStart = "ATG"

// IsRNAStart is an effecient case insensitve test for the RNA start codon.
func IsRNAStart(b0, b1, b2 byte) bool {
	// force bytes to upper case
	return b0&0xdf == RNAStart[0] &&
		b1&0xdf == RNAStart[1] &&
		b2&0xdf == RNAStart[2]
}

// IsDNAStart is an effecient case insensitve test for the DNA start codon.
func IsDNAStart(b0, b1, b2 byte) bool {
	// force bytes to upper case
	return b0&0xdf == DNAStart[0] &&
		b1&0xdf == DNAStart[1] &&
		b2&0xdf == DNAStart[2]
}

var dnaStartRx = regexp.MustCompile("[Aa][Tt][Gg]")

// Translate translates RNA and returns an amino acid sequence.
//
// Translation begins at the start of the sequence.  The sequence should begin
// with the RNA start codon.  Translation ends at a stop codon or at the
// end of the sequence.
//
// Errors are returned for no start codon and no stop codon but the
// translated sequence is returned in any case.
func (s RNA8) Translate() (a AA20, err error) {
	if len(s) < 3 || !IsRNAStart(s[0], s[1], s[2]) {
		err = errors.New("No start codon")
	}
	return translate8(s)
}

func translate8(s []byte) (a AA20, err error) {
	for p := 0; ; p += 3 {
		if p+3 > len(s) {
			err = errors.New("No stop codon")
			return
		}
		aa := translateCodon8(s[p], s[p+1], s[p+2])
		if aa == aaStop {
			return
		}
		a = append(a, aa)
	}
}

// TranslateORF locates and translates all open reading frames in a sequence.
//
// Returned is a collection of all unique amino acid sequences.
func (s DNA8) TranslateORF() (c []AA20) {
	m := map[string]AA20{}
	orf := func(s DNA8) {
		for {
			start := dnaStartRx.FindIndex(s)
			if start == nil {
				return
			}
			s = s[start[0]:]
			ps, err := translate8(s)
			if err == nil {
				m[string(ps)] = ps
			}
			s = s[1:]
		}
	}
	orf(s)
	orf(s.ReverseComplement())
	for _, s := range m {
		c = append(c, s)
	}
	return
}

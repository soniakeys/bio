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

// CodonTable allows lookup of an amino acid symbol by a 6 bit index
// corresponding to a codon.
//
// Use CodonIndex() to compute a suitable index.
//
// Use the function TranslateCodon though in preference to directly indexing
// the table.
var CodonTable [64]byte

func init() {
	copy(CodonTable[:], []byte(
		"KNNKTTTTIIIMRSSRQHHQPPPPLLLLRRRR*YY*SSSSLFFL*CCWEDDEAAAAVVVVGGGG"))
}

const AAStop = '*' // symbol for stop codon

// CodonIndex computes a 6 bit index suitable for indexing CodonTable.
//
// Arguments are three consecutive symbols of RNA or DNA and should be
// DNA8 or RNA8 symbols.
func CodonIndex(b0, b1, b2 byte) int {
	// &6 bit mask produces these values:
	// A   = 00
	// C   = 01
	// U,T = 10
	// G   = 11
	return int(b0&6<<3 | b1&6<<1 | b2&6>>1)
}

// TranslateCodon translates a sequence of three DNA8 or RNA8 bases into
// an amino acid symbol of AA20Alphabet or the stop symbol AAStop.
func TranslateCodon(b0, b1, b2 byte) byte {
	return CodonTable[CodonIndex(b0, b1, b2)]
}

const (
	DNAStart = "ATG" // DNA start codon
	RNAStart = "AUG" // RNA start codon
)

// IsStartCodon tests for a DNA or RNA start codon.
//
// The arguments should be a sequence of DNA8 or RNA8 symbols.
func IsStartCodon(b0, b1, b2 byte) bool {
	return b0&6 == 'a'&6 && b1&6 == 'u'&6 && b2&6 == 'g'&6
}

// Translate8 translates DNA or RNA into an amino acid sequence.
//
// Argument seq should be DNA8 or RNA8.  Sequences containing other symbols
// will give nonsense results.
//
// Translation begins at the beginning of the sequence.  A start codon
// is not required.  Translation ends at a stop codon or at the end of
// the sequence.
//
// An error is returned if there is no stop codon.
// The translated sequence is always returned.
func Translate8(seq []byte) (AA20, error) {
	var a AA20
	for p := 2; p < len(seq); p += 3 {
		aa := TranslateCodon(seq[p-2], seq[p-1], seq[p])
		if aa == AAStop {
			return a, nil
		}
		a = append(a, aa)
	}
	return a, errors.New("No stop codon")
}

/*
// translate8 translates a DNA8 or RNA8 sequence into an AA20 sequence.
//
// The function can translate either DNA or RNA.  The sequence s must begin with
// the specified start codon which should be either the DNAStart or RNAStart.
// Translation ends at a stop codon, which must be present.
// Errors are returned for no start codon or no stop codon.
func translate8(s, startCodon []byte) (AA20, error) {
	if !bytes.HasPrefix(s, startCodon) {
		return nil, errors.New("No start codon")
	}
	a := AA20{'M'}
	for p := 3; ; p += 3 {
		if p+3 > len(s) {
			return nil, errors.New("No stop codon")
		}
		aa := CodonTable[CodonIndex(s[p], s[p+1], s[p+2])]
		if aa == AAStop {
			break
		}
		a = append(a, aa)
	}
	return a, nil
}
*/
var dnaStartRx = regexp.MustCompile("[Aa][Tt][Gg]")

// TranslateORF locates and translates all open reading frames in a sequence.
//
// Returned is a collection of all unique amino acid sequences.
func (s DNA8) TranslateORF() []AA20 {
	// TODO binary insertion with sort.Search would be more efficient than
	// a map here.
	m := map[string]AA20{}
	orf := func(s DNA8) {
		for {
			start := dnaStartRx.FindIndex(s)
			if start == nil {
				return
			}
			s = s[start[0]:]
			ps, err := Translate8(s)
			if err == nil {
				m[string(ps)] = ps
			}
			s = s[1:]
		}
	}
	orf(s)
	orf(s.ReverseComplement())
	c := make([]AA20, len(m))
	i := 0
	for _, s := range m {
		c[i] = s
		i++
	}
	return c
}

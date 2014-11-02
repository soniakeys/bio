// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

// translate.go
//
// Definitions concerning tranlation of DNA and RNA into amino acid sequences.

package bio

import (
	"bytes"
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
var codonSlice = CodonTable[:]
var codonBMap = map[[3]byte]byte{}
var codonSMap = map[string]byte{}
var codonInv ['Z'][][]byte   // codons indexed by amino acid
var codonInvRC ['Z'][][]byte // RC of codons indexed by amino acid
var codonInvRx ['Z']string
var codonInvRCRx ['Z']string

func init() {
	copy(CodonTable[:], []byte(
		"KNNKTTTTIIIMRSSRQHHQPPPPLLLLRRRR*YY*SSSSLFFL*CCWEDDEAAAAVVVVGGGG"))
	c := DNA8("AAA")
	var a [3]byte
	for _, aa := range CodonTable {
		copy(a[:], c)
		codonBMap[a] = aa
		codonSMap[string(c)] = aa
		if aa != AAStop {
			codonInv[aa] = append(codonInv[aa], append([]byte{}, c...))
			codonInvRC[aa] = append(codonInvRC[aa], c.ReverseComplement())
		}
		c.Inc()
	}
	for _, aa := range AA20Alphabet {
		codonInvRx[aa] =
			"((" + string(bytes.Join(codonInv[aa], []byte(")|("))) + "))"
		codonInvRCRx[aa] =
			"((" + string(bytes.Join(codonInvRC[aa], []byte(")|("))) + "))"
	}
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

func TranslateCodonC(b0, b1, b2 byte) byte {
	return codonSlice[CodonIndex(b0, b1, b2)]
}

func TranslateCodonB(b [3]byte) byte {
	return codonBMap[b]
}

func TranslateCodonS(s string) byte {
	return codonSMap[s]
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

// translate8 translates DNA or RNA into an amino acid sequence.
//
// Argument seq should be DNA8 or RNA8.  Sequences containing other symbols
// will give nonsense results.
func translate8(seq []byte) (a AA20, stop bool) {
	for p := 2; p < len(seq); p += 3 {
		aa := TranslateCodon(seq[p-2], seq[p-1], seq[p])
		if aa == AAStop {
			return a, true
		}
		a = append(a, aa)
	}
	return a, false
}

// Translate translates the receiver DNA sequence into an amino acid
// sequence.
//
// Translation begins at the beginning of the sequence.  A start codon
// is not required.  Translation ends at a stop codon or at the end of
// the sequence.
//
// The translated sequence is always returned regardless of the presence
// of a stop codon.  Return value stop is true if a stop codon was present.
func (s DNA8) Translate() (p AA20, stop bool) { return translate8(s) }

// Translate translates the receiver RNA sequence into an amino acid
// sequence.
//
// Translation begins at the beginning of the sequence.  A start codon
// is not required.  Translation ends at a stop codon or at the end of
// the sequence.
//
// The translated sequence is always returned regardless of the presence
// of a stop codon.  Return value stop is true if a stop codon was present.
func (s RNA8) Translate() (p AA20, stop bool) { return translate8(s) }

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
			ps, stop := s.Translate()
			if stop {
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

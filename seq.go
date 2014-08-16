// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"bytes"
	"errors"
)

// seq.go
//
// Functions general to multiple sequence types.
//
// There is no named type for generalized sequences.  []byte seems enough.
// Functions here are generally case sensitive and should be documented as
// such, to avoid confusion with the behavior of case insensitive types
// such as DNA8.

const LCBit = 0x20 // Bit mask for ASCII alphabetic values

// Freq returns counts of all symbols appearing in a sequence.
//
// The returned map contains a key for each symbol appearing in the sequence.
// If a symbol does not appear in the sequence, there will be no map key for it.
// Upper and lower case are treated as distinct.  If a sequence contains both
// 'A' and 'a' the result will contain separate counts for the two symbols.
func Freq(s []byte) map[byte]int {
	m := map[byte]int{}
	for _, b := range s {
		m[b]++
	}
	return m
}

// case insensitive.  used by methods of both DNA8 and RNA8.
func baseFreq8(s []byte) (a, c, tu, g int) {
	var n [8]int
	for _, b := range s {
		n[b&6]++
	}
	return n['A'&6], n['C'&6], n['T'&6], n['G'&6]
}

// Hamming returns the Hamming distance between two byte sequences.
// Comparison is done byte-wise and so is case sensitive.
func Hamming(s, t []byte) (int, error) {
	if len(t) != len(s) {
		return 0, errors.New("Hamming: unequal lengths")
	}
	h := 0
	for i, b := range s {
		if b != t[i] {
			h++
		}
	}
	return h, nil
}

// AllIndex finds all occurrences of a motif in a sequence.
//
// Returned is a list of indexes of all occurrences of motif m in sequence s,
// including overlapping ones.   Comparison is done byte-wise and so is
// case sensitive.
func AllIndex(s, m []byte) (x []int) {
	for searched := 0; ; {
		i := bytes.Index(s[searched:], m)
		if i < 0 {
			break
		}
		searched += i + 1
		x = append(x, searched)
	}
	return
}

// ToLower returns a new sequence with ASCII upper case forced to lower case.
//
// Non-alphabetic bytes are unaffected.
func ToLower(s []byte) []byte {
	r := make([]byte, len(s))
	for i, b := range s {
		if b >= 'A' && b <= 'Z' {
			b |= LCBit
		}
		r[i] = b
	}
	return r
}

// ToUpper returns a new sequence with ASCII lower case forced to upper case.
//
// Non-alphabetic bytes are unaffected.
func ToUpper(s []byte) []byte {
	r := make([]byte, len(s))
	for i, b := range s {
		if b >= 'a' && b <= 'z' {
			b &^= LCBit
		}
		r[i] = b
	}
	return r
}

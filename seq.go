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
// such as a reminder that are different than methods on base and amino acid
// types.

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

// function common to both DNA8 and RNA8.
func baseFreq8(s []byte) (a, c, tu, g int) {
	var n [4]int
	for _, b := range s {
		n[b>>1&3]++
	}
	return n['A'>>1&3], n['C'>>1&3], n['T'>>1&3], n['G'>>1&3]
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

// ToLower returns a new sequence with upper case symbols forced to lower case.
//
// This is the wrong way to convert case in natural languge text that could
// include non-ASCII letters, but it is the right way to convert sequences of
// byte-length symbols, like base sequences for example.  Use this function
// as needed to prepare sequences for case sensitive functions.
func ToLower(s []byte) []byte {
	r := make([]byte, len(s))
	for i, b := range s {
		if b >= 'A' && b <= 'Z' {
			b += 32
		}
		r[i] = b
	}
	return r
}

// ToUpper returns a new sequence with lower case symbols forced to upper case.
//
// This is the wrong way to convert case in natural languge text that could
// include non-ASCII letters, but it is the right way to convert sequences of
// byte-length symbols, like base sequences for example.  Use this function
// as needed to prepare sequences for case sensitive functions.
func ToUpper(s []byte) []byte {
	r := make([]byte, len(s))
	for i, b := range s {
		if b >= 'a' && b <= 'z' {
			b -= 32
		}
		r[i] = b
	}
	return r
}

// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

// Package dna defines operations on DNA strings.
//
// It defines types and methods that are specific to DNA or optimized for DNA.
package dna

// dna.go
//
// Types Seq, Strict, and related definitions.

// Seq type represents a DNA string.
//
// It is expected to hold DNA base symbols but other symbols are allowed.
// Methods on the type accomodate non-base symbols.
type Seq []byte

// Strict type represents a string consisting strictly of DNA symbols.
//
// Allowed symbols are ACTGactg.  Methods on the type assume this.  They
// are case-insensitive but may produce nonsense results if the string
// contains non-base symbols.
type Strict []byte

// String satisfies fmt.Stringer.
func (s Seq) String() string {
	return string(s)
}

// String satisfies fmt.Stringer.
func (s Strict) String() string {
	return string(s)
}

// Freq returns counts of all symbols appearing in the string.
//
// The returned map contains a key for each symbol appearing in the string.
// If a symbol does not appear in the string, there will be no map key for it.
// Upper and lower case are treated as distinct.  If a string contains both
// 'A' and 'a' the result will contain separate counts for the two symbols.
func (s Seq) Freq() map[byte]int {
	m := map[byte]int{}
	for _, b := range s {
		m[b]++
	}
	return m
}

// BaseFreq returns counts of each of the four DNA bases.
func (s Strict) BaseFreq() (a, c, t, g int) {
	var n [4]int
	for _, b := range s {
		n[b>>1&3]++
	}
	return n['A'>>1&3], n['C'>>1&3], n['T'>>1&3], n['G'>>1&3]
}

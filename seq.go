// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

// seq.go
//
// Functions general to multiple sequence types.
//
// There is no named type for generalized sequences.  []byte seems enough.

// Freq returns counts of all symbols appearing in a string.
//
// The returned map contains a key for each symbol appearing in the string.
// If a symbol does not appear in the string, there will be no map key for it.
// Upper and lower case are treated as distinct.  If a string contains both
// 'A' and 'a' the result will contain separate counts for the two symbols.
func Freq(s []byte) map[byte]int {
	m := map[byte]int{}
	for _, b := range s {
		m[b]++
	}
	return m
}

// function common to both DNAStrict and RNAStrict.
func baseFreq(s []byte) (a, c, tu, g int) {
	var n [4]int
	for _, b := range s {
		n[b>>1&3]++
	}
	return n['A'>>1&3], n['C'>>1&3], n['T'>>1&3], n['G'>>1&3]
}

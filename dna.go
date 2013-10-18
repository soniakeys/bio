// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

// dna.go
//
// Types and methods that are specific to DNA or optimized for DNA.

// DNA type represents a DNA string.
//
// It is expected to hold DNA base symbols but other symbols are allowed.
// Methods on the type accomodate non-base symbols.
type DNA []byte

// DNAStrict type represents a string consisting strictly of DNA symbols.
//
// Allowed symbols are ACTGactg.  Methods on the type assume this.  Methods
// are case-insensitive but may produce nonsense results if the string
// contains non-base symbols.
type DNAStrict []byte

// String satisfies fmt.Stringer.
func (s DNA) String() string {
	return string(s)
}

// String satisfies fmt.Stringer.
func (s DNAStrict) String() string {
	return string(s)
}

// BaseFreq returns counts of each of the four DNA bases.
func (s DNAStrict) BaseFreq() (a, c, t, g int) {
	return baseFreq(s)
}

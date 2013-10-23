// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
//	"errors"
)

// rna.go
//
// Types and methods that are specific to RNA or optimized for RNA.

// RNA type represents a RNA sequence.
//
// It is expected to hold RNA base symbols but other symbols are allowed.
// Methods on the type accomodate non-base symbols and accomodate both
// upper and lower case base symbols.
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

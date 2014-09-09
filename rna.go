// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

// rna.go
//
// Types and methods that are specific to RNA or optimized for RNA.

// RNA type represents a RNA sequence.
//
// It is expected to hold RNA base symbols but other symbols are allowed.
// Methods on the type accomodate non-base symbols and accomodate both
// upper and lower case base symbols.
type RNA []byte

// String converts the receiver to a string.
func (s RNA) String() string {
	return string(s)
}

// RNA8 type represents a sequence of upper or lower case RNA symbols.
//
// Allowed symbols are the eight symbols ACUGacug.  Methods on the type
// assume this.  Methods are thus case-insensitive but may produce nonsense
// results if the sequence contains non-base symbols.
type RNA8 []byte

// String satisfies fmt.Stringer.
func (s RNA8) String() string {
	return string(s)
}

// BaseFreq returns counts of each of the four RNA bases.
func (s RNA8) BaseFreq() (a, c, u, g int) {
	return baseFreq8(s)
}

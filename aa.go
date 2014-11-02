// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

// aa.go
//
// Amino acid definitions

// AA type holds amino acid sequences.
//
// Symbols are not restricted to AA20Alphabet, data may include other
// symbols such a gap symbols.
type AA []byte

// String returns an AA converted to a string.
func (s AA) String() string {
	return string(s)
}

// AA20 type holds amino acid sequences.  Content should be strictly limited
// to the symbols of AA20Alphabet, in upper case.  Methods may panic
// on other symbols.
type AA20 []byte

var AA20Alphabet = AA20("ACDEFGHIKLMNPQRSTVWY") // IUPAC symbols for the 20 proteinogenic amino acids

// String returns an AA20 converted to a string.
func (s AA20) String() string {
	return string(s)
}

// NumSubPepLinear computes the theoretical number of subpeptides of a
// linear peptide of length l.
func NumSubPepLinear(l int) int {
	return (l*(l+1))/2 + 1
}

// NumSubPepCyclic computes the theoretical number of subpeptides of a
// cyclic peptide of length l.
func NumSubPepCyclic(l int) int {
	return l * (l - 1)
}

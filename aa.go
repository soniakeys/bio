// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

// aa.go
//
// Amino acid definitions

// AA type holds amino acid sequences.
type AA []byte

// AA20 holds amino acid sequences.  The content is strictly limited
// to the twenty proteinogenic symbols, in upper case.
type AA20 []byte

// AA20Alphabet is the set of 20 proteinogenic amino acid symbols.
const AA20Alphabet = "ACDEFGHIKLMNPQRSTVWY"

func (s AA) String() string {
	return string(s)
}

func (s AA20) String() string {
	return string(s)
}

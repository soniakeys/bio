// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

// aa.go
//
// Amino acid definitions

// AA type holds amino acid sequences.
type AA []byte

// AAStrict holds amino acid sequences.  The content is strictly limited
// to the 20 proteinogenic symbols.
type AAStrict []byte

// AAAlphabet is the set of 20 proteinogenic amino acid symbols.
const AAAlphabet = "ACDEFGHIKLMNPQRSTVWY"

func (s AA) String() string {
	return string(s)
}

func (s AAStrict) String() string {
	return string(s)
}

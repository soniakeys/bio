// Bio is a package of bioinformatics algorithms.
//
// Bio is an experimental project for exploring concepts in bioinformatics.
// It is not currently intended for any practical use.
//
// Data types for sequences
//
// The concept is to use very simple types with the goal of performance,
// testing the idea that lots of generality might not be needed.
//
// Byte slices, as []byte with no named type are used as the most general
// sequnce type and a number of functions will operate on this type.  These
// functions will generally implement some string algorithm that is not
// specific to DNA or peptides.  These functions will generally treat []byte
// as a sequence of symbols and perform case sensitive symbol comparisons.
//
// A number of types are defined as named types for byte slices and are
// intended to hold data restricted to some symbol set.  These are named
// with some indication of the expected symbol set and are documented with
// a list or description of the symbols.  Methods on the type will assume
// the symbol set.  For example the type DNA8 holds the eight symbols
// ACGTacgt.
//
// The package API does not duplicate methods across all possible types.
// The main purpose of the package is to explore algorithms, not to be
// comprehensive.  Further, it's not obvious that orthogonality would even
// be necessary in practical applications.
package bio

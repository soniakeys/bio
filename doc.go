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
// ACTGactg.
//
// The package API does not duplicate methods across all possible types.
// The main purpose of the package is to explore algorithms, not to be
// comprehensive.  Further, it's not obvious that orthogonality would even
// be necessary in practical applications.
//
// Bit twiddling and base order
//
// This library usually works with bases in the order ACTG and ACUG.
// In this order, bits 1 and 2 of the ASCII representations are 00, 01, 10,
// and 11, allowing indexes 0, 1, 2, 3 to be directly obtained from the ASCII.
// Yes, T and U are the same, 10.  In some cases this bit pattern is exploited
// for efficiency, in other cases this base order is maintained just for
// consistency within the library.
//
// Reference implementations
//
// Much of the code in this library is experimental or derived from coursework
// exercises.  Some though is implemented from pseudocode in published sources.
// This code will be commented with one of the following references:
//
// Jones 2004: An Introduction to Bioinformatics Algorithms, Neil C. Jones
// and Pavel A. Pevzner, MIT Press, 2004.
//
// Rajasekaran 2005: "Exact Algorithms for Planted Motif Problems",
// Sanguthevar Rajasekaran, Journal of Computational Biology, 2005
//
// Compeau 2014: Bioinformatics Algorithms, An Active Learning Approach,
// Phillip Compeau and Pavel Pevzner, Active Learning Publishers, 2014.
package bio

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

func transcribe(s []byte) []byte {
	t := append([]byte{}, s...)
	for i, b := range t {
		if b&0xdf == 'T' {
			t[i]++
		}
	}
	return t
}

// Transcribe returns the RNA transcription of the receiver.
//
// A new string is returned.  The receiver string is unmodified.
func (s DNA) Transcribe() RNA {
	return RNA(transcribe(s))
}

// Transcribe returns the RNA transcription of the receiver.
//
// A new string is returned.  The original string is unmodified.
func (s DNAStrict) Transcribe() RNAStrict {
	return RNAStrict(transcribe(s))
}

// ReverseComplement returns the reverse complement of the receiver.
//
// A new string is returned.  The receiver is left unmodified.
// Symbols not in the DNA alphabet are reversed but otherwise left unchanged.
func (s DNA) ReverseComplement() DNA {
	rc := make(DNA, len(s))
	rcx := len(rc)
	for _, b := range s {
		switch b | 0x20 {
		case 'a', 't':
			b ^= 0x15
		case 'c', 'g':
			b ^= 0x04
		}
		rcx--
		rc[rcx] = b
	}
	return rc
}

// ReverseComplement returns the reverse complement of the receiver.
//
// A new string is returned.  The receiver is left unmodified.
func (s DNAStrict) ReverseComplement() DNAStrict {
	rc := make(DNAStrict, len(s))
	rcx := len(rc)
	for _, b := range s {
		rcx--
		rc[rcx] = ^b&2>>1*17 | 4 ^ b
		// it was faster, at least on the machine it was written on.
	}
	return rc
}

// CGFraction returns the fraction of the string that is C or G over the
// string length.
func (s DNAStrict) CGFraction() float64 {
	_, c, _, g := baseFreq(s)
	return float64(c+g) / float64(len(s))
}

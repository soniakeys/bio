// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

// dna.go
//
// Types and methods that are specific to DNA or optimized for DNA.

// DNA type represents a DNA sequence.
//
// DNA is expected to generally hold DNA base symbols but other symbols are
// allowed.  Methods on the type process both upper and lower case symbols
// for the four DNA bases.  Methods accomodate other symbols but generally
// do not process them except as documented.
type DNA []byte

// DNA8 type represents a sequence of upper or lower case DNA symbols.
//
// Allowed symbols are the eight symbols ACTGactg.  Methods on the type
// assume this.  Methods are thus case-insensitive but may produce nonsense
// results if the sequence contains other symbols.
type DNA8 []byte

// String satisfies fmt.Stringer.
func (s DNA) String() string {
	return string(s)
}

// String satisfies fmt.Stringer.
func (s DNA8) String() string {
	return string(s)
}

// BaseFreq returns counts of each of the four DNA bases.
func (s DNA8) BaseFreq() (a, c, t, g int) {
	return baseFreq8(s)
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
// A new sequence is returned.  The receiver is unmodified.
func (s DNA) Transcribe() RNA {
	return RNA(transcribe(s))
}

// Transcribe returns the RNA transcription of the receiver.
//
// A new sequence is returned.  The receiver is unmodified.
func (s DNA8) Transcribe() RNA8 {
	return RNA8(transcribe(s))
}

// ReverseComplement returns the reverse complement of the receiver.
//
// A new sequence is returned.  The receiver is unmodified.
// Case is maintainted for symbols in the DNA alphabet.
// Symbols not in the DNA alphabet are reversed but otherwise left unchanged.
func (s DNA) ReverseComplement() DNA {
	rc := make(DNA, len(s))
	rcx := len(rc)
	for _, b := range s {
		switch b | LCBit {
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
// A new sequence is returned.  The receiver is unmodified.
func (s DNA8) ReverseComplement() DNA8 {
	rc := make(DNA8, len(s))
	rcx := len(rc)
	for _, b := range s {
		rcx--
		rc[rcx] = ^b&2>>1*17 | 4 ^ b
		// crazy bit twiddling, but it was faster,
		// at least on the machine it was written on.
	}
	return rc
}

// GCContent returns the fraction of the sequence that is C or G over the
// string length.  The value returned is in the range 0 to 1.
func (s DNA8) GCContent() float64 {
	_, c, _, g := baseFreq8(s)
	return float64(c+g) / float64(len(s))
}

// DNAConsensus returns a consensus sequence from multiple sequences.
//
// Consensus in each position is simply the most frequent base in that
// position.  If, for a given position, a base does not appear in any sequence,
// a '-' is emitted.  Input sequences should be of the same lengths, but the
// result will be the length of the first sequence.  Sequences shorter or
// longer than the first are allowed, any excess length being ignored.  While
// the function is case insensitive, the result is returned as upper case.
//
// Score is the sum of occurrences of the consensus base at each position
// over the sequence.  Maximum possible score is len(c) * len(c[0]), which
// would happen if all sequences were identical.
func DNAConsensus(c []DNA) (seq DNA, score int) {
	if len(c) == 0 {
		return
	}
	s := c[0]
	if len(s) == 0 {
		return
	}
	const bases = "A C T G"
	r := make(DNA, len(s))
	// profile posistion by position, without constructing profile matrix
	for i := range r {
		// profile
		var n [7]int
		for _, s := range c {
			if i < len(s) {
				// c.f. DNA8 version
				switch b := s[i] | LCBit; b {
				case 'a', 'c', 't', 'g':
					n[b&6]++
				}

			}
		}
		// find consensus
		max := n[0]
		maxb := 0
		for b := 2; b < 8; b += 2 {
			if n[b] > max {
				max = n[b]
				maxb = b
			}
		}
		// c.f. DNA8 version
		if max > 0 {
			r[i] = bases[maxb]
			score += max
		} else {
			r[i] = '-'
		}
	}
	return r, score
}

// DNA8Consensus returns a consensus sequence from multiple sequences.
//
// Consensus in each position is simply the most frequent base in that
// position.  Input sequences should be of the same lengths, but the result
// will be the length of the first sequence.  Sequences shorter or longer
// than the first are allowed, any excess length being ignored.  While the
// function is case insensitive, the result is returned as upper case.
//
// Score is the sum of occurrences of the consensus base at each position
// over the sequence.  Maximum possible score is len(c) * len(c[0]), which
// would happen if all sequences were identical.
func DNA8Consensus(c []DNA8) (seq DNA8, score int) {
	if len(c) == 0 {
		return
	}
	s := c[0]
	if len(s) == 0 {
		return
	}
	const bases = "A C T G"
	r := make(DNA8, len(s))
	// profile posistion by position, without constructing profile matrix
	for i := range r {
		// profile
		var n [7]int
		for _, s := range c {
			if i < len(s) {
				// c.f. DNA version
				n[s[i]&6]++
			}
		}
		// find consensus
		max := n[0]
		maxb := 0
		for b := 2; b < 8; b += 2 {
			if n[b] > max {
				max = n[b]
				maxb = b
			}
		}
		// c.f. DNA version
		r[i] = bases[maxb]
		score += max
	}
	return r, score
}

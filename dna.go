// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

// dna.go
//
// Types and methods that are specific to DNA or optimized for DNA.

// DNA type represents a DNA sequence.
//
// It is expected to hold DNA base symbols but other symbols are allowed.
// Methods on the type accomodate non-base symbols and accomodate
// both upper and lower case base symbols.
type DNA []byte

// DNAStrict type represents a sequence consisting strictly of DNA symbols.
//
// Allowed symbols are ACTGactg.  Methods on the type assume this.  Methods
// are case-insensitive but may produce nonsense results if the sequence
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
// A new sequence is returned.  The receiver is unmodified.
func (s DNA) Transcribe() RNA {
	return RNA(transcribe(s))
}

// Transcribe returns the RNA transcription of the receiver.
//
// A new sequence is returned.  The receiver is unmodified.
func (s DNAStrict) Transcribe() RNAStrict {
	return RNAStrict(transcribe(s))
}

// ReverseComplement returns the reverse complement of the receiver.
//
// A new sequence is returned.  The receiver is unmodified.
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
// A new sequence is returned.  The receiver is unmodified.
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

// CGFraction returns the fraction of the sequence that is C or G over the
// string length.
func (s DNAStrict) CGFraction() float64 {
	_, c, _, g := baseFreq(s)
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
	const bases = "ACTG"
	r := make(DNA, len(s))
	// profile posistion by position, without constructing profile matrix
	for i := range r {
		// profile
		var n [4]int
		for _, s := range c {
			if i < len(s) {
				// c.f. DNAStrict version
				switch b := s[i] | 0x20; b {
				case 'a', 'c', 't', 'g':
					n[b>>1&3]++
				}

			}
		}
		// find consensus
		max := n[0]
		maxb := 0
		for b := 1; b < 4; b++ {
			if n[b] > max {
				max = n[b]
				maxb = b
			}
		}
		// c.f. DNAStrict version
		if max > 0 {
			r[i] = bases[maxb]
			score += max
		} else {
			r[i] = '-'
		}
	}
	return r, score
}

// DNAStrictConsensus returns a consensus sequence from multiple sequences.
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
func DNAStrictConsensus(c []DNAStrict) (seq DNAStrict, score int) {
	if len(c) == 0 {
		return
	}
	s := c[0]
	if len(s) == 0 {
		return
	}
	const bases = "ACTG"
	r := make(DNAStrict, len(s))
	// profile posistion by position, without constructing profile matrix
	for i := range r {
		// profile
		var n [4]int
		for _, s := range c {
			if i < len(s) {
				// c.f. DNA version
				n[s[i]>>1&3]++
			}
		}
		// find consensus
		max := n[0]
		maxb := 0
		for b := 1; b < 4; b++ {
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

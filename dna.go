// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

// DNA type represents a DNA sequence.
//
// DNA is expected to generally hold DNA base symbols but other symbols are
// allowed.  Methods on the type process both upper and lower case symbols
// for the four DNA bases.  Methods accommodate other symbols but generally
// do not process them except as documented.
type DNA []byte

// String satisfies fmt.Stringer.
func (s DNA) String() string {
	return string(s)
}

// BaseFreq returns the counts of each of the four DNA bases.
//
// Symbols which are not DNA bases are ignored and not included in any count.
func (s DNA) BaseFreq() (a, c, t, g int) {
	f := Seq(s).Freq()
	return f['A'] + f['a'],
		f['C'] + f['c'],
		f['T'] + f['t'],
		f['G'] + f['g']
}

// Transcribe returns the RNA transcription of the receiver.
//
// A new sequence is returned.  The receiver is unmodified.
func (s DNA) Transcribe() RNA {
	t := append(RNA{}, s...)
	for i, b := range t {
		if b|LCBit == 't' {
			t[i]++
		}
	}
	return t
}

var iupacDNAComp [256]byte

func init() {
	iupacDNAComp = [256]byte{
		'G': 'C',
		'T': 'A',
		'A': 'T',
		'C': 'G',
		'R': 'Y',
		'Y': 'R',
		'S': 'S',
		'W': 'W',
		'K': 'M',
		'M': 'K',
		'D': 'H',
		'H': 'D',
		'B': 'V',
		'V': 'B',
		'N': 'N',
		'g': 'c',
		't': 'a',
		'a': 't',
		'c': 'g',
		'r': 'y',
		'y': 'r',
		's': 's',
		'w': 'w',
		'k': 'm',
		'm': 'k',
		'd': 'h',
		'h': 'd',
		'b': 'v',
		'v': 'b',
		'n': 'n',
	}
	for i, b := range iupacDNAComp {
		if b == 0 {
			iupacDNAComp[i] = byte(i)
		}
	}
}

// DNAComplement returns the complement of a DNA symbol.
//
// DNAComplement complements all IUPAC DNA symbols, including ambiguity symbols.
// It allows lower case versions IUPAC symbols and returns the complement
// with case preserved.  If the symbol is not an IUPAC DNA symbol,
// it is returned unchanged.
func DNAComplement(b byte) byte {
	return iupacDNAComp[b]
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
		rcx--
		rc[rcx] = DNAComplement(b)
	}
	return rc
}

// GCContent returns the fraction of the sequence that is G or C over
// all DNA bases in the sequence, ignoring case.
// The value returned is in the range 0 to 1.
func (s DNA) GCContent() float64 {
	a, c, t, g := s.BaseFreq()
	gc := float64(c + g)
	return gc / (gc + float64(a+t))
}

// TiTvRatio computes the transition to transversion ratio of two
// DNA strings.
//
// Non-DNA symbols are ignored.
//
// The strings must be of equal length.  Panic or nonsense result
// if the strings are of unequal length.
func TiTvRatio(s, t DNA) float64 {
	var ti, tv int
	for i, s1 := range s {
		s1 |= LCBit
		switch t1 := t[i] | LCBit; t1 {
		case s1:
			continue
		case 'a', 'g':
			switch s1 {
			case 'a', 'g':
				ti++
			case 'c', 't':
				tv++
			}
		case 'c', 't':
			switch s1 {
			case 'c', 't':
				ti++
			case 'a', 'g':
				tv++
			}
		}
	}
	return float64(ti) / float64(tv)
}

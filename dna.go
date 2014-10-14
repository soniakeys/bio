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

// BaseFreq returns the counts of each of the four DNA bases.
//
// Symbols which are not DNA bases are ignored and not included in any count.
func (s DNA) BaseFreq() (a, c, t, g int) {
	f := Freq(s)
	return f['A'] + f['a'],
		f['C'] + f['c'],
		f['T'] + f['t'],
		f['G'] + f['g']
}

// BaseFreq returns counts of each of the four DNA bases.
func (s DNA8) BaseFreq() (a, c, t, g int) {
	return baseFreq8(s)
}

func transcribe(s []byte) []byte {
	t := append([]byte{}, s...)
	for i, b := range t {
		if b|LCBit == 't' {
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

// DNAComplement returns the complement of a DNA symbol.  It complements
// DNA base symbols, preserving case.  If the symbol is not a DNA symbol,
// it is returned unchanged.
func DNAComplement(b byte) byte {
	switch b | LCBit {
	case 'a', 't':
		return b ^ 0x15
	case 'c', 'g':
		return b ^ 0x04
	}
	return b
}

// DNA8Complement returns the complement of a DNA8 symbol preserving case.
// If the symbol is not a DNA8 symbol, the result is nonsense.
func DNA8Complement(b byte) byte {
	return ^b&2>>1*17 | 4 ^ b
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

// ReverseComplement returns the reverse complement of the receiver.
//
// A new sequence is returned.  The receiver is unmodified.
func (s DNA8) ReverseComplement() DNA8 {
	rc := make(DNA8, len(s))
	rcx := len(rc)
	for _, b := range s {
		rcx--
		rc[rcx] = DNA8Complement(b)
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

// GCContent returns the fraction of the sequence that is G or C over
// the string length.  The value returned is in the range 0 to 1.
func (s DNA8) GCContent() float64 {
	_, c, _, g := baseFreq8(s)
	return float64(c+g) / float64(len(s))
}

// GCSkew senses G and C in a DNA8 symbol.
//
// Result:
//   1 for 'G' or 'g'
//   -1 for 'C' or 'c'
//   0 for 'A', 'a', 'T' or 't'
//   Nonsense for other symbols
func GCSkew(b byte) int {
	b >>= 1
	return int(int8(-(b & 1) & (b&2 - 1)))
}

// MinGCSkew returns indexes in s with the minimum cumulative GC skew
// from the beginning of the string.
func (s DNA8) MinGCSkew() (m []int) {
	min := 0
	skew := 0
	for i, b := range s {
		skew += GCSkew(b)
		switch {
		case skew < min:
			m = []int{i}
			min = skew
		case skew == min:
			m = append(m, i)
		}
	}
	return
}

// HammingVariants returns a list of all DNA k-mers within hamming distance h
// of receiver kmer k.  Case is preserved by position.
func (k DNA8) HammingVariants(h int) []DNA8 {
	v := []DNA8{append(DNA8{}, k...)}
	const sym = "A C T G"
	var f func(DNA8, int)
	f = func(t DNA8, h int) {
		for i := 0; i < len(t); i++ {
			sub := t[i:]
			b := sub[0]
			vb := byte(0)
			for j := 0; j < 3; j++ {
				if vb == b&6 {
					vb += 2
				}
				sub[0] = sym[vb] | b&32
				v = append(v, append(DNA8{}, k...))
				if h > 1 && len(sub) > 1 {
					f(sub[1:], h-1)
				}
				vb += 2
			}
			sub[0] = b
		}
	}
	f(k, h)
	return v
}

// ModalHammingKmers returns DNA kmers matching a maximum number of times
// to string s, where matching kmers match some kmer of s within hamming
// distance d.
func (s DNA8) ModalHammingKmers(k, d int) (m []DNA8) {
	c := map[string][]DNA8{}
	f := map[string]int{}
	max := 0
	for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
		k0 := s[i:j]
		s0 := string(k0)
		v, ok := c[s0]
		if !ok {
			v = k0.HammingVariants(d)
			c[s0] = v
		}
		for _, kmer := range v {
			s := string(kmer)
			n := f[s] + 1
			f[s] = n
			switch {
			case n > max:
				m = []DNA8{kmer}
				max = n
			case n == max:
				m = append(m, kmer)
			}
		}
	}
	return
}

// ModalHammingKmersRC returns DNA kmers matching a maximum number of times
// to string s, where matching kmers or their reverse complements match some
// kmer of s within hamming distance d.
func (s DNA8) ModalHammingKmersRC(k, h int) (m []DNA8) {
	c := map[string][]DNA8{}
	f := map[string]int{}
	max := 0
	// TODO pull this out and share with ModalHammingKmers
	tally := func(k0 DNA8) {
		s0 := string(k0)
		v, ok := c[s0]
		if !ok {
			v = k0.HammingVariants(h)
			c[s0] = v
		}
		for _, kmer := range v {
			s := string(kmer)
			n := f[s] + 1
			f[s] = n
			switch {
			case n > max:
				m = []DNA8{kmer}
				max = n
			case n == max:
				m = append(m, kmer)
			}
		}
	}
	r := s.ReverseComplement()
	for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
		tally(s[i:j])
		tally(r[i:j])
	}
	return
}

// Find indexes in s where s translates to pep.  Searches all three
// reading frames and finds overlaps but does not search reverse complement.
func (s DNA8) AAFindAllIndex(pep AA20) (r []int) {
	t := make(AA20, len(s)/3)
	for f := 0; f < 3; f++ {
		t = t[:(len(s)-f)/3]
		for i, j := f, 0; j < len(t); i, j = i+3, j+1 {
			t[j] = TranslateCodon(s[i], s[i+1], s[i+2])
		}
		x := AllIndex(t, pep)
		for i, p := range x {
			x[i] = f + p*3
		}
		r = append(r, x...)
	}
	return
}

// Find indexes in s where s or reverse complement of s translates to pep.
// Searches all six reading frames, finds overlaps.  Returns 0-based indexes
// from the start of s.
func (s DNA8) AAFindAllIndexRC(pep AA20) []int {
	f := s.AAFindAllIndex(pep)
	r := s.ReverseComplement().AAFindAllIndex(pep)
	for i, p := range r {
		r[i] = len(s) - p - len(pep)*3
	}
	return append(f, r...)
}

// Hamming returns hamming distance between s and t.  Nonsense or panic
// results if strings are unequal length.
func (s DNA8) Hamming(t DNA8) (d int) {
	for i, b := range s {
		if b&6 != t[i]&6 {
			d++
		}
	}
	return
}

// MotifSeqDist returns the minimum hamming distance from motif m
// to any same length kmer in sequence s.
func (m DNA8) MotifSeqDist(s DNA8) int {
	min := len(m)
	for i, j := 0, len(m); j < len(s); i, j = i+1, j+1 {
		if h := m.Hamming(s[i:j]); h < min {
			min = h
		}
	}
	return min
}

// MotifSetDist is a distance measure from a motif m to a set
// of strings l.
//
// (Not a mathematical set, just a list.)
//
// It is the sum of distances MotifSeqDist from m to each string in l.
func (m DNA8) MotifSetDist(l []DNA8) int {
	d := 0
	for _, s := range l {
		d += m.MotifSeqDist(s)
	}
	return d
}

// KmersNearestMotif returns the kmers in s having minimum hamming
// distance from motif m.
func (m DNA8) KmersNearestMotif(s DNA8) (k []DNA8) {
	min := len(m)
	for i, j := 0, len(m); j < len(s); i, j = i+1, j+1 {
		switch h := m.Hamming(s[i:j]); {
		case h < min:
			min = h
			k = []DNA8{s[i:j]}
		case h == min:
			k = append(k, s[i:j])
		}
	}
	return
}

// Inc "increments" a kmer, for the purpose of iterating over all possible
// kmers.  The symbol order is ACTG.  A string of all Gs rolls over to
// all As.
func (m DNA8) Inc() {
	for i := len(m) - 1; i >= 0; i-- {
		b := m[i]
		if n := b & 6; n < 6 {
			m[i] = "C T G"[n] | b&LCBit
			return
		}
		m[i] = 'A' | b&LCBit
	}
}

// TTRatio compultes the transition to transversion ratio of two
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

// TTRatio compultes the transition to transversion ratio of two
// DNA8 strings.
//
// The strings must be of equal length.  Panic or nonsense result
// if the strings are of unequal length or if non-DNA8 symbols are present.
func TiTvRatio8(s, t DNA8) float64 {
	var ti, tv int
	for i, s1 := range s {
		s1 &= 6
		switch t1 := t[i] & 6; {
		case s1+t1 == 6:
			ti++
		case s1 != t1:
			tv++
		}
	}
	return float64(ti) / float64(tv)
}

// PalIndex is a return type for FindPalAllIndex.
type PalIndex struct {
	Index int // index where a palindrome was found
	Len   int // length of the palindrome
}

// PalFindAllIndex finds palindrome sequences.
//
// A palindrome sequence equals its reverse complement.  Arguments min and max
// are minimum and maximum lengths of palindrome sequences to find.
//
// Returned is a list of indexes and lengths of all palindrome sequences
// with length from min to max inclusive.
func (s DNA8) PalFindAllIndex(min, max int) (p []PalIndex) {
	switch {
	case min < 2:
		min = 2
	case min&1 == 1:
		min++
	}
	if len(s) < min {
		return nil
	}
	if len(s) < max {
		max = len(s)
	}
	if max&1 == 1 {
		max--
	}
	if max < min {
		return nil
	}

	// the algorithm scans a sliding window of the receiver sequence s.
	// the window size starts at min for the beginning of the seqence and
	// increments to max for scanning the interior of the sequence, then
	// decrements back down to min at the end.  actually the window size
	// isn't directly needed -- more useful is inner, the index of the inner
	// pair in t.
	inner := min/2 - 1 // initial value corresponding to min window size.

	x := 0 // index of window in s

	// scan from inner and work out.
	// if a check fails at any point, skip remaining checks.
	// for all palindromes >= min, append result.
	scan := func() {
		for i, j, n := inner, inner+1, 2; i >= 0; i, j, n = i-1, j+1, n+2 {
			// bit twiddling solution by "i"
			xr := (s[x+i] ^ s[x+j]) &^ LCBit
			if xr != 'A'^'T' && xr != 'C'^'G' {
				break
			}
			if n >= min {
				p = append(p, PalIndex{x + i, n})
			}
		}
	}
	// short windows at the beginning of s
	for imax := max/2 - 1; inner < imax; inner++ {
		scan()
	}
	// interior of s
	for end := len(s) - max; x < end; x++ {
		scan()
	}
	// windows at end of s
	for end := len(s) - min; x <= end; x += 2 {
		scan()
		inner--
	}
	return p
}

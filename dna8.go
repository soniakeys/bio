// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"bytes"
	"math/big"
	"regexp"
)

// DNA8 type represents a sequence of upper or lower case DNA symbols.
//
// Allowed symbols are the eight symbols ACTGactg.  Methods on the type
// assume this.  Methods are thus case-insensitive but may produce nonsense
// results if the sequence contains other symbols.
type DNA8 []byte

// String satisfies fmt.Stringer.
func (s DNA8) String() string {
	return string(s)
}

// Cmp is a comparison function, suitable for sorting.
//
// Rules:
//
// 1. A shorter sequence is < a longer sequence.
//
// 2. The order of bases is ACTG.
//
// 3. A lower case base is < the upper case of the same base.
//
// 4. As usual for the DNA8 type, results for non-DNA8 data are undefined.
//
// Returns -1 for s < t, 0 for s == t, 1 for s > t
func (s DNA8) Cmp(t DNA8) int {
	// rule 1:
	switch {
	case len(s) < len(t):
		return -1
	case len(s) > len(t):
		return 1
	}
	for i, si := range s {
		ti := t[i]
		switch {
		// rule 2:
		case si&6 < ti&6:
			return -1
		case si&6 > ti&6:
			return 1
		// rule 3:
		case si < ti:
			return 1
		case si > ti:
			return -1
		}
	}
	return 0
}

// AllIndex finds all occurrences of a motif in a sequence.
//
// Returned is a list of indexes of all occurrences of motif m in sequence s,
// including overlapping ones.
func (s DNA8) AllIndex(m DNA8) (x []int) {
ij:
	for i, j := 0, len(m); j <= len(s); i, j = i+1, j+1 {
		for k, mk := range m {
			if (s[i+k]^mk)&^LCBit != 0 {
				continue ij
			}
		}
		x = append(x, i)
	}
	return
}

// AllCount counts all occurrences of a motif in a sequence,
// including overlaps.
//
// It is equivalent to len(s.AllIndex(m)).
func (s DNA8) AllCount(m DNA8) (x int) {
ij:
	for i, j := 0, len(m); j <= len(s); i, j = i+1, j+1 {
		for k, mk := range m {
			if (s[i+k]^mk)&^LCBit != 0 {
				continue ij
			}
		}
		x++
	}
	return
}

// ModalKmers returns the most frequent k-mers in a DNA8 string.
//
//  s is the string to search.
//  k is the k-mer length.
//
// This variant works for large k, but at the expense of allocated copy of s.
// See ModalSmallKmers for a variant more efficent with small k.
func (s DNA8) ModalKmers(k int) Kmers {
	ms := Str(s).ModalKmers(k) // call string version; it uses a map.
	m := make(Kmers, len(ms))
	for i, mi := range ms {
		m[i] = DNA8(mi)
	}
	return m
}

func (s DNA8) num() int64 {
	var o int64
	for _, b := range s {
		o = o<<2 | int64(b&6)
	}
	return o >> 1
}

func strDNA8(n int64, k int) DNA8 {
	a := make(DNA8, k)
	for k > 0 {
		k--
		a[k] = "ACTG"[n&3]
		n >>= 2
	}
	return a
}

// FreqArray constructs a "frequency array" representing frequencies of kmers
// occurring in s.
//
// The returned slice will have length 4^k representing frequencies of kmers
// in lexicographic order where the alphabet order is ACTG.
//
// Compare to function KmerComposition which also computes kmer frequencies.
// Freq array is somewhat more efficient but only practical with small k.
func (s DNA8) FreqArray(k int) []int {
	a := make([]int, 1<<(2*uint(k)))
	mask := int64(len(a) - 1)
	n := s[:k].num()
	a[n] = 1
	for k < len(s) {
		n = n<<2&mask | int64(s[k]>>1&3)
		a[n]++
		k++
	}
	return a
}

// ModalSmallKmers returns the most frequent k-mers in a DNA8 string.
//
//  s is the string to search.
//  k is the k-mer length.
//
// This variant is fast and efficient--as long as k isn't too big.
// It allocates no additional copy of s, but does allocate an integer array
// of size 2^(2*k).  So it's good with large s and small k.
// See ModalKmers for a variant practical for larger k.
func (s DNA8) ModalSmallKmers(k int) Kmers {
	a := s.FreqArray(k)
	var max, nMax int
	for _, f := range a {
		switch {
		case f == max:
			nMax++
		case f > max:
			max = f
			nMax = 1
		}
	}
	m := make(Kmers, 0, nMax)
	for n, f := range a {
		if f == max {
			m = append(m, strDNA8(int64(n), k))
		}
	}
	return m
}

// BaseFreq returns counts of each of the four DNA bases.
func (s DNA8) BaseFreq() (a, c, t, g int) {
	return baseFreq8(s)
}

// Transcribe returns the RNA transcription of the receiver.
//
// A new sequence is returned.  The receiver is unmodified.
func (s DNA8) Transcribe() RNA8 {
	return RNA8(DNA(s).Transcribe())
}

// DNA8Complement returns the complement of a DNA8 symbol preserving case.
// If the symbol is not a DNA8 symbol, the result is nonsense.
func DNA8Complement(b byte) byte {
	return ^b&2>>1*17 | 4 ^ b
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

// NumHammingVariants computes the expected length of the result of
// DNA8.HammingVariants.
func NumHammingVariants(k, d int) *big.Int {
	if d > k {
		d = k
	}
	k64 := int64(k)
	d64 := int64(d)
	three := big.NewInt(3)
	p3 := big.NewInt(1) // a power of 3
	sum := big.NewInt(1)
	var b big.Int
	for h := int64(1); h <= d64; h++ {
		b.Binomial(k64, h)
		sum.Add(sum, b.Mul(&b, p3.Mul(p3, three)))
	}
	return sum
}

// HammingVariants returns a list of all DNA k-mers within hamming distance d
// of receiver kmer k.  Case is preserved by position.
func (kmer DNA8) HammingVariants(d int) Kmers {
	// recursive, but minimizes allocations
	if d == 0 {
		return Kmers{append(DNA8{}, kmer...)}
	}
	v := Kmers{append(DNA8{}, kmer...)}
	const sym = "A C T G"
	var f func(DNA8, int)
	f = func(t DNA8, h int) {
		for i := 0; i < len(t); i++ {
			sub := t[i:]
			b := sub[0]
			for vb := byte(0); ; vb += 2 {
				if vb == b&6 {
					vb += 2
				}
				if vb == 8 {
					break
				}
				sub[0] = sym[vb] | b&32
				v = append(v, append(DNA8{}, kmer...))
				if h > 1 && len(sub) > 1 {
					f(sub[1:], h-1)
				}
			}
			sub[0] = b
		}
	}
	f(kmer, d)
	return v
}

// HammingVariants1.  Same result, different algorithm.
func (kmer DNA8) HammingVariantsRef(d int) Kmers {
	// "by the book."  well, except it's improved to preserve case by position.
	// churns memory a bit.
	if d == 0 {
		return Kmers{append(DNA8{}, kmer...)}
	}
	if len(kmer) == 1 {
		b := kmer[0]
		c := b & LCBit
		return Kmers{
			DNA8{'A' | c},
			DNA8{'C' | c},
			DNA8{'T' | c},
			DNA8{'G' | c},
		}
	}
	var nd Kmers
	suf := kmer[1:]
	c := kmer[0] & LCBit
	for _, txt := range suf.HammingVariants(d) {
		if suf.Hamming(txt) < d {
			for _, b := range DNA8("ACTG") {
				nd = append(nd, append(DNA8{b | c}, txt...))
			}
		} else {
			nd = append(nd, append(DNA8{kmer[0]}, txt...))
		}
	}
	return nd
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

// Four algorithms for AAFindAllIndex.  Typical benchmarks:
//
// bio $ go test -run AAFind -bench AAFind
// PASS
// BenchmarkAAFindAllIndex-8      	     100	  28360381 ns/op
// BenchmarkAAFindAllIndex3-8     	      50	  27133149 ns/op
// BenchmarkAAFindAllIndexOnline-8	      30	  37121451 ns/op
// BenchmarkAAFindAllIndexRegexp-8	       2	 682268192 ns/op
//
// Benchmarks in dna8_test.go but usually commented out.

// AAFindAllIndex returnes indexes in s where a substring of s or the
// reverse complement of the substring translates to pep.
//
// It searches all six reading frames and includes overlaps and returns a
// list of 0-based indexes from the start of s.  Unfortunately results are
// not necessarily ordered and translation palindromes are included twice.
//
// Compared to the alternative algorithms this one is both simple and gives
// good performance.
func (s DNA8) AAFindAllIndex(pep AA20) []int {
	f := s.aaFindAllIndex(pep)
	r := s.ReverseComplement().aaFindAllIndex(pep)
	for i, p := range r {
		r[i] = len(s) - p - len(pep)*3
	}
	return append(f, r...)
}

// Helper method for AAFindAllIndex.  Searches all three
// reading frames and finds overlaps but does not search reverse complement.
func (s DNA8) aaFindAllIndex(pep AA20) (r []int) {
	t := make(AA20, len(s)/3)
	for f := 0; f < 3; f++ {
		t = t[:(len(s)-f)/3]
		for i, j := f, 0; j < len(t); i, j = i+3, j+1 {
			t[j] = TranslateCodon(s[i], s[i+1], s[i+2])
		}
		x := Seq(t).AllIndex(Seq(pep))
		for i, p := range x {
			x[i] = f + p*3
		}
		r = append(r, x...)
	}
	return
}

// AAFindAllIndex3 is an alternative algorithm for AAFindAllIndex.
//
// Like AAFindAllIndex, results are not necessarily in order and translation
// palindromes are counted twice.
//
// The algorithm is somewhat more complex and runs only slightly faster
// than AAFindAllIndex.
func (s DNA8) AAFindAllIndex3(pep AA20) []int {
	f := s.aaFindAllIndex3(pep)
	r := s.ReverseComplement().aaFindAllIndex3(pep)
	for i, p := range r {
		r[i] = len(s) - p - len(pep)*3
	}
	return append(f, r...)
}

// helper function for AAFindAllIndex3
func (s DNA8) aaFindAllIndex3(pep AA20) (r []int) {
	if len(s) < 3 {
		return
	}
	lt := (len(s) + 2) / 3
	t0 := make(AA20, lt*3)
	t2 := t0[2*lt:]
	t1 := t0[lt : 2*lt]
	t0 = t0[:lt]
	cx := s[0]&6<<1 | s[1]&6>>1
	i := 2
	for tx := 0; ; tx++ {
		cx = cx<<2 | s[i]&6>>1
		t0[tx] = CodonTable[cx&63]
		i++
		if i == len(s) {
			break
		}

		cx = cx<<2 | s[i]&6>>1
		t1[tx] = CodonTable[cx&63]
		i++
		if i == len(s) {
			break
		}

		cx = cx<<2 | s[i]&6>>1
		t2[tx] = CodonTable[cx&63]
		i++
		if i == len(s) {
			break
		}
	}
	x0 := Seq(t0).AllIndex(Seq(pep))
	for i, p := range x0 {
		x0[i] = p * 3
	}
	x1 := Seq(t1).AllIndex(Seq(pep))
	for i, p := range x1 {
		x1[i] = p*3 + 1
	}
	x2 := Seq(t2).AllIndex(Seq(pep))
	for i, p := range x2 {
		x2[i] = p*3 + 2
	}
	return append(x0, append(x1, x2...)...)
}

// AAFindAllIndexOnline is an alternative algorithm for AAFindAllIndex.
//
// It is an online algorithm that makes a single pass over s and uses
// additional memory proportional to |pep| rather than |s|.
//
// It has the additional advantage that the result list of indexes comes
// out in order, and without repeating indexes of translation palindromes.
//
// Sadly the added complexity of the algorithm seems to leave it a little
// (maybe 25%) slower than AAFindAllIndex.
func (s DNA8) AAFindAllIndexOnline(pep AA20) (x []int) {
	// algorithm: use six ring buffers, three for forward reading frames and
	// three for reverse.  handle each codon just once, in order.  for each
	// codon, translate to a single position in a forward ring buffer and
	// reverse translate to a single position in a reverse ring buffer.
	// fill all but one position in each of the six buffers as a priming step,
	// then begin update-test cycle.  store the match if either forward or
	// reverse buffers match the peptide.
	var f, r [3][]byte // ring buffers for the six reading frames
	for i := range f {
		f[i] = make([]byte, len(pep))
		r[i] = make([]byte, len(pep))
	}
	p3 := 3 * len(pep) // size of window into Text
	if p3 > len(s) {
		return // not enough dna
	}
	// prime buffers.  i is start of codon being translated.
	prEnd := p3 - 3
	for i := 0; i < prEnd; i++ {
		fx := i % 3             // reading frame index
		px := i / 3             // peptide position index
		rx := len(pep) - 1 - px // revc peptide index
		f[fx][px] = TranslateCodon(s[i], s[i+1], s[i+2])
		r[fx][rx] = TranslateCodonC(s[i+2], s[i+1], s[i])
	}
	// scan remainder of s, testing for peptide
	sEnd := len(s) - 3
	for i := prEnd; i <= sEnd; i++ {
		fx := i % 3             // reading frame index
		px := i / 3 % len(pep)  // peptide position index
		rx := len(pep) - 1 - px // revc peptide index
		f[fx][px] = TranslateCodon(s[i], s[i+1], s[i+2])
		r[fx][rx] = TranslateCodonC(s[i+2], s[i+1], s[i])
		fMatch := bytes.Equal(f[fx][px+1:], pep[:rx]) &&
			bytes.Equal(f[fx][:px+1], pep[rx:])
		rMatch := bytes.Equal(r[fx][:rx], pep[px+1:]) &&
			bytes.Equal(r[fx][rx:], pep[:px+1])
		if fMatch || rMatch {
			x = append(x, i-prEnd)
		}
	}
	return
}

// AAFindAllIndexRegexp is an alternative algorithm to AAFindAllIndex.
//
// It is slower.
func (s DNA8) AAFindAllIndexRegexp(pep AA20) (r []int) {
	var pat, rcPat string
	for _, aa := range pep {
		pat += codonInvRx[aa]
		rcPat = codonInvRCRx[aa] + rcPat
	}
	pat = "(" + pat + ")|(" + rcPat + ")"
	rx := regexp.MustCompile(pat)
	for searched := 0; ; {
		p := rx.FindIndex(s[searched:])
		if p == nil {
			break
		}
		x := searched + p[0]
		r = append(r, x)
		searched = x + 1
	}
	return r
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

/*
// MotifHamming returns the minimum hamming distance from motif m
// to any same length kmer in sequence s.
func (s DNA8) MotifHamming(m DNA8) int {
	min := len(m)
	for i, j := 0, len(m); j <= len(s); i, j = i+1, j+1 {
		if h := m.Hamming(s[i:j]); h < min {
			min = h
		}
	}
	return min
}
*/
// MotifHamming returns the minimum hamming distance from motif m
// to any same length kmer in sequence s.
func (s DNA8) MotifHamming(m DNA8) int {
	min := len(m)
pos:
	for i, j := 0, len(m); j <= len(s); i, j = i+1, j+1 {
		h := 0
		for i, b := range s[i:j] {
			if m[i] != b {
				if h == min {
					continue pos
				}
				h++
			}
		}
		min = h
	}
	return min
}

// MotifHamming is a distance measure from a motif m to a list of sequences l.
//
// It is the sum of distances DNA8.MotifHamming for string in l.
//
// Sequences in receiver l may be of different lengths.
func (l DNA8List) MotifHamming(m DNA8) int {
	d := 0
	for _, s := range l {
		d += s.MotifHamming(m)
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
// kmers.  The symbol order is ACTG.  Sequences < all Gs increment and return
// true.  A sequence of all Gs rolls over to all As and returns false.
func (m DNA8) Inc() bool {
	for i := len(m) - 1; i >= 0; i-- {
		b := m[i]
		if n := b & 6; n < 6 {
			m[i] = "C T G"[n] | b&LCBit
			return true
		}
		m[i] = 'A' | b&LCBit
	}
	return false
}

// nextVertex and bypass are used by DNAList.MedianMotifsB.  They seem a
// little quirky and hard to explain so they're not exported for now.
func (a DNA8) nextVertex(i int) int {
	if i < len(a) {
		a[i] = 'A' | a[i]&LCBit
		return i + 1
	}
	for j := len(a) - 1; j >= 0; j-- {
		b := a[j]
		if n := b & 6; n < 6 {
			a[j] = "C T G"[n] | b&LCBit
			return j + 1
		}
	}
	return 0
}

func (a DNA8) bypass(i int) int {
	for j := i - 1; j >= 0; j-- {
		b := a[j]
		if n := b & 6; n < 6 {
			a[j] = "C T G"[n] | b&LCBit
			return j + 1
		}
	}
	return 0
}

// TiTvRatio8 compultes the transition to transversion ratio of two
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

func (s DNA8) UniqueKmers(k int) map[Str]struct{} {
	t := Str(Seq(s).ToUpper())
	m := map[Str]struct{}{}
	for i, j := 0, k; j <= len(t); i, j = i+1, j+1 {
		m[t[i:j]] = struct{}{}
	}
	return m
}

func (s DNA8) UniqueHammingKmers(k, d int) map[Str]struct{} {
	h := map[Str]struct{}{}
	var vs StrKmers
	for u := range s.UniqueKmers(k) {
		vs = u.DNA8HammingVariants(d, vs[:0])
		for _, v := range vs {
			h[v] = struct{}{}
		}
	}
	return h
}

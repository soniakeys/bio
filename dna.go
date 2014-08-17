// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"bytes"
	"math"
	"math/rand"
	"sort"
)

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
func (s DNA) BaseFreq() (a, c, g, t int) {
	f := Freq(s)
	return f['A'] + f['a'],
		f['C'] + f['c'],
		f['G'] + f['g'],
		f['T'] + f['t']
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
	a, c, g, t := s.BaseFreq()
	gc := float64(c + g)
	return gc / (gc + float64(a+t))
}

// GCContent returns the fraction of the sequence that is G or C over
// the string length.  The value returned is in the range 0 to 1.
func (s DNA8) GCContent() float64 {
	_, c, _, g := baseFreq8(s)
	return float64(c+g) / float64(len(s))
}

// DNAProfileMatrix represents a DNA profile matrix, useful for constructing
// consensus sequences.
//
// Construct with make(DNAProfileMatrix, seqLen)
type DNAProfileMatrix []struct{ A, C, G, T int }

// Add adds DNA sequence s to the profile matrix.  Mixed case is allowed.
// String lengths not equal to the profile matrix length are allowed,
// with shortages and excesses ignored.
// Symbols other than DNA bases are ignored.
func (pm DNAProfileMatrix) Add(s DNA) {
	if len(s) > len(pm) {
		s = s[:len(pm)]
	}
	for x, b := range s {
		switch b | LCBit {
		case 'a':
			pm[x].A++
		case 'c':
			pm[x].C++
		case 'g':
			pm[x].G++
		case 't':
			pm[x].T++
		}
	}
}

// Consensus returns a consensus string from a populated DNA profile matrix.
// If no valid DNA symbol ocurred in a particular position in all input
// sequences, a '-' is emitted in that position.
func (pm DNAProfileMatrix) Consensus() DNA {
	cs := make(DNA, len(pm))
	for x := range cs {
		c := byte('-') // consensus symbol
		fc := 0        // freq of consensus symbol
		if f := pm[x].A; f > fc {
			c = 'A'
			fc = f
		}
		if f := pm[x].C; f > fc {
			c = 'C'
			fc = f
		}
		if f := pm[x].G; f > fc {
			c = 'G'
			fc = f
		}
		if pm[x].T > fc {
			c = 'T'
		}
		cs[x] = c
	}
	return cs
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
				// (see DNA8 version)
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
		// (see DNA8 version)
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
				// (see DNA version)
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
		// (see DNA version)
		r[i] = bases[maxb]
		score += max
	}
	return r, score
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

// MinGCSkew returns the positions in s with the minimum cumulative GC skew
// from the beginning of the string.
func (s DNA8) MinSkew() (m []int) {
	min := 0
	skew := 0
	for i, b := range s {
		skew += GCSkew(b)
		switch {
		case skew < min:
			m = []int{i + 1}
			min = skew
		case skew == min:
			m = append(m, i+1)
		}
	}
	return
}

// HammingVariants returns a list of all DNA k-mers within hamming distance h
// of receiver kmer k.  Case is preserved by position.
func (k DNA8) HammingVariants(h int) []string {
	v := []string{string(k)}
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
				v = append(v, string(k))
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

// ModalVariantKmer returns the most frequent DNA k-mers within hamming
// distance h of any k-mer present in receiver string s.
func (s DNA8) ModalVariantKmer(k, h int) (m []string) {
	c := map[string][]string{}
	f := map[string]int{}
	max := 0
	for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
		k0 := s[i:j]
		s0 := string(k0)
		v, ok := c[s0]
		if !ok {
			v = k0.HammingVariants(h)
			c[s0] = v
		}
		for _, kmer := range v {
			n := f[kmer] + 1
			f[kmer] = n
			switch {
			case n > max:
				m = []string{kmer}
				max = n
			case n == max:
				m = append(m, kmer)
			}
		}
	}
	return
}

// ModalVariantKmer returns the most frequent DNA k-mers within hamming
// distance h of any k-mer present in either receiver string s or the
// reverse complement of s
func (s DNA8) ModalVariantKmerRC(k, h int) (m []string) {
	c := map[string][]string{}
	f := map[string]int{}
	max := 0
	// TODO pull this out and share with ModalVariantKmer
	tally := func(k0 DNA8) {
		s0 := string(k0)
		v, ok := c[s0]
		if !ok {
			v = k0.HammingVariants(h)
			c[s0] = v
		}
		for _, kmer := range v {
			n := f[kmer] + 1
			f[kmer] = n
			switch {
			case n > max:
				m = []string{kmer}
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
		x := AllIndex(pep, t)
		for i, p := range x {
			x[i] = p*3 + f
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

// MotifMatrix is a list of kmers all with the same length.
//
// Methods may return nonsense results or panic if kmers are not of the same
// length.
type MotifMatrix []DNA8

// Score returns the number of unconserved bases in matrix x.  A base is
// conserved when it is the modal base in a position.  If multimodal, only
// a single base is considered conserved.
func (x MotifMatrix) Score() (s int) {
	k := len(x[0])
	for i := 0; i < k; i++ {
		var f [4]int
		max := 0
		for _, m := range x {
			b := m[i] >> 1 & 3
			f[b]++
			if f[b] > max {
				max = f[b]
			}
		}
		s += len(x) - max
	}
	return
}

// Count returns a 4xk matrix of the counts of each base.
// The row order is ACTG.
func (x MotifMatrix) Count() (actg [4][]int) {
	k := len(x[0])
	f := [4][]int{
		make([]int, k),
		make([]int, k),
		make([]int, k),
		make([]int, k)}
	for _, m := range x {
		for i, b := range m {
			f[b>>1&3][i]++
		}
	}
	return f
}

// Profile represents a profile matrix corresponding to a MotifMatrix.
// Elements are probabilities and sum to 1 in each column.
// Row order is ACTG.
type Profile [4][]float64

// NewProfile allocates a profile matrix, leaving it all zeros.
func NewProfile(k int) Profile {
	return Profile{
		make([]float64, k),
		make([]float64, k),
		make([]float64, k),
		make([]float64, k)}
}

// NewPseudoProfile allocates a profile matrix, initializing it with
// pseudocount probabilty p.
func NewPseudoProfile(k int, p float64) Profile {
	a := make([]float64, k*4)
	for i := range a {
		a[i] = p
	}
	return Profile{
		a[:k],
		a[k : k*2],
		a[k*2 : k*3],
		a[k*3:]}
}

// KmerProb computes the probability of a kmer given profile p.
func (p Profile) KmerProb(kmer DNA8) float64 {
	pr := 1.
	for i, b := range kmer {
		pr *= p[b>>1&3][i]
	}
	return pr
}

// Kmer returns the profile-most probable kmers in s.
func (p Profile) Kmer(s DNA8) (k []DNA8) {
	max := 0.
	for i, j := 0, len(p[0]); j < len(s); i, j = i+1, j+1 {
		switch pr := p.KmerProb(s[i:j]); {
		case pr > max:
			k = []DNA8{s[i:j]}
			max = pr
		case pr == max:
			k = append(k, s[i:j])
		}
	}
	return
}

func (p Profile) RandKmer(s DNA8) DNA8 {
	k := len(p[0])
	c := make([]float64, len(s)-k+1)
	c[0] = p.KmerProb(s[:k])
	for i, j := 1, k+1; j < len(s); i, j = i+1, j+1 {
		c[i] = c[i-1] + p.KmerProb(s[i:j])
	}
	i := sort.SearchFloat64s(c, c[len(c)-1]*rand.Float64())
	return s[i : i+k]
}

func GibbsSampler(l []DNA8, k, N int) (MotifMatrix, int) {
	motifs := make(MotifMatrix, len(l))
	for i, s := range l {
		j := rand.Intn(len(s) - k + 1)
		motifs[i] = s[j : j+k]
	}
	best := len(l) * k
	bestMotifs := make(MotifMatrix, len(l))
	for i := 0; i < N; i++ {
		j := rand.Intn(len(l))
		motifs[j] = motifs[0]
		motifs[0], motifs[j] = motifs[j], motifs[1:].Profile().RandKmer(l[j])
		if s := motifs.Score(); s < best {
			best = s
			copy(bestMotifs, motifs)
		}
	}
	return bestMotifs, best
}

func (p Profile) Motifs(l []DNA8) MotifMatrix {
	m := make(MotifMatrix, len(l))
	for i, s := range l {
		m[i] = p.Kmer(s)[0]
	}
	return m
}

func RandomMotifSearch(l []DNA8, k int) (MotifMatrix, int) {
	motifs := make(MotifMatrix, len(l))
	for i, s := range l {
		j := rand.Intn(len(s) - k + 1)
		motifs[i] = s[j : j+k]
	}
	bestMotifs := make(MotifMatrix, len(l))
	best := len(l) * k
	for {
		motifs = motifs.LaplaceProfile().Motifs(l)
		if s := motifs.Score(); s < best {
			copy(bestMotifs, motifs)
			best = s
		} else {
			return bestMotifs, s
		}
	}
}

// Profile constructs the profile matrix corresponding to motif matrix m.
func (x MotifMatrix) Profile() Profile {
	p := NewProfile(len(x[0]))
	inc := 1 / float64(len(x))
	for _, m := range x {
		for i, b := range m {
			p[b>>1&3][i] += inc
		}
	}
	return p
}

// LaplaceProfile constructs the profile matrix corresponding to
// motif matrix m, augmented by Laplace's Rule of Succession.  That is,
// with a pseudocount of 1 added to each base count.
func (x MotifMatrix) LaplaceProfile() Profile {
	inc := 1 / float64(len(x)+4)
	p := NewPseudoProfile(len(x[0]), inc)
	for _, m := range x {
		for i, b := range m {
			p[b>>1&3][i] += inc
		}
	}
	return p
}

// Consensus generates a consensus string from motif matrix x.
func (x MotifMatrix) Consensus() DNA8 {
	k := len(x[0])
	c := make(DNA8, k)
	for i := 0; i < k; i++ {
		var f [4]int
		maxf := 0
		var maxb byte
		for _, m := range x {
			b := m[i] >> 1 & 3
			f[b]++
			if f[b] > maxf {
				maxf = f[b]
				maxb = b
			}
		}
		c[i] = "ACTG"[maxb]
	}
	return c
}

// Entropy computes entropy for a motif matrix.  It is the sum of entropies
// in each column.
func (x MotifMatrix) Entropy() (e float64) {
	k := len(x[0])
	inc := 1 / float64(len(x))
	for i := 0; i < k; i++ {
		var f [4]float64
		for _, m := range x {
			f[m[i]>>1&3] += inc
		}
		for _, p := range f {
			if p > 0 {
				e -= p * math.Log2(p)
			}
		}
	}
	return
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
func (m DNA8) MotifSeqDistance(s DNA8) int {
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
		d += m.MotifSeqDistance(s)
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
	for i, b := range m {
		if n := b & 6; n < 6 {
			m[i] = "C T G"[n] | b&32
			return
		}
		m[i] = 'A' | b&32
	}
}

// MedianString returns a list of kmers that are at minimum distance
// (by the method MotifSetDist) to a list of strings l.
// The algorithm is brute force and practical only when k is small.
func MedianString(l []DNA8, k int) (m []DNA8) {
	z := DNA8(bytes.Repeat([]byte{'A'}, k))
	min := len(l[0])
	for p := append(DNA8{}, z...); ; {
		switch d := p.MotifSetDist(l); {
		case d < min:
			m = []DNA8{append(DNA8{}, p...)}
			min = d
		case d == min:
			m = append(m, append(DNA8{}, p...))
		default:
		}
		p.Inc()
		if bytes.Equal(p, z) {
			break
		}
	}
	return
}

// MotifSearch uses a greedy algorithm.
func MotifSearch(l []DNA8, k int) (m []DNA8) {
	bestMotifs := make(MotifMatrix, len(l))
	for i, s := range l {
		bestMotifs[i] = s[:k]
	}
	motifs := make(MotifMatrix, len(l))
	s0 := l[0]
	best := len(motifs) * k
	for i, j := 0, k; j < len(s0); i, j = i+1, j+1 {
		motifs[0] = s0[i:j]
		for i := 1; i < len(l); i++ {
			p := motifs[:i].LaplaceProfile()
			motifs[i] = p.Kmer(l[i])[0]
		}
		if s := motifs.Score(); s < best {
			best = s
			copy(bestMotifs, motifs)
		}
	}
	return bestMotifs
}

// TTRatio compultes the translation to transversion ratio of two
// DNA strings.
//
// The strings must be of equal length.
//
// Non-DNA symbols are ignored, but the function returns false
// if the strings are of unequal length.
func TTRatio(s, t DNA) (float64, bool) {
	if len(t) != len(s) {
		return 0, false
	}
	var ts, tv int
	for i, si := range s {
		switch ti := t[i]; ti {
		case si:
			continue
		case 'A', 'G':
			switch si {
			case 'A', 'G':
				ts++
			case 'C', 'T':
				tv++
			}
		case 'C', 'T':
			switch si {
			case 'C', 'T':
				ts++
			case 'A', 'G':
				tv++
			}
		}
	}
	return float64(ts) / float64(tv), true
}

// EqRev4 and EqRevSpan are a bit quirky.
//
// EqRev4 tests if the first four symbols of s are the reverse
// of the first four symbols of t, ignoring case
//
// It panics if the sequences do not have at least four symbols.
func EqRev4(s, t DNA) bool {
	return EqRevIndex(s, t, 3) && EqRevIndex(s[1:], t[1:], 1)
}

// EqRevIndex returns true if symbols in two sequences are equal
// but reversed.
//
// It tests the symbols s[0] against t[index] and test s[index] against t[0].
// The comparison is case insensitive.
// The function panics if the strings do not have at least index+1 symbols.
func EqRevIndex(s, t DNA, index int) bool {
	return (s[0]^t[index]|s[index]^t[0])&^LCBit == 0
}

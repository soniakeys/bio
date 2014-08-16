// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"bytes"
	"sort"
	"strings"
)

// seq.go
//
// Functions general to multiple sequence types.
//
// There is no named type for generalized sequences.  []byte seems enough.
// Functions here are generally case sensitive and should be documented as
// such, to avoid confusion with the behavior of case insensitive types
// such as DNA8.

const LCBit = 0x20 // Bit mask for ASCII alphabetic values

// Freq returns counts of all symbols appearing in a sequence.
//
// The returned map contains a key for each symbol appearing in the sequence.
// If a symbol does not appear in the sequence, there will be no map key for it.
// Upper and lower case are treated as distinct.  If a sequence contains both
// 'A' and 'a' the result will contain separate counts for the two symbols.
func Freq(s []byte) map[byte]int {
	m := map[byte]int{}
	for _, b := range s {
		m[b]++
	}
	return m
}

// case insensitive.  used by methods of both DNA8 and RNA8.
func baseFreq8(s []byte) (a, c, tu, g int) {
	var n [8]int
	for _, b := range s {
		n[b&6]++
	}
	return n['A'&6], n['C'&6], n['T'&6], n['G'&6]
}

// Hamming returns the Hamming distance between two byte sequences.
// Comparison is done byte-wise and so is case sensitive.
//
// Sequences s and t must be of equal length.
// Panic or nonsense results if the lengths are unequal.
func Hamming(s, t []byte) int {
	h := 0
	for i, b := range s {
		if b != t[i] {
			h++
		}
	}
	return h
}

// AllIndex finds all occurrences of a motif in a sequence.
//
// Returned is a list of indexes of all occurrences of motif m in sequence s,
// including overlapping ones.   Comparison is done byte-wise and so is
// case sensitive.
func AllIndex(s, m []byte) (x []int) {
	for searched := 0; ; {
		i := bytes.Index(s[searched:], m)
		if i < 0 {
			break
		}
		searched += i + 1
		x = append(x, searched)
	}
	return
}

// ToLower returns a new sequence with ASCII upper case forced to lower case.
//
// Non-alphabetic bytes are unaffected.
func ToLower(s []byte) []byte {
	r := make([]byte, len(s))
	for i, b := range s {
		if b >= 'A' && b <= 'Z' {
			b |= LCBit
		}
		r[i] = b
	}
	return r
}

// ToUpper returns a new sequence with ASCII lower case forced to upper case.
//
// Non-alphabetic bytes are unaffected.
func ToUpper(s []byte) []byte {
	r := make([]byte, len(s))
	for i, b := range s {
		if b >= 'a' && b <= 'z' {
			b &^= LCBit
		}
		r[i] = b
	}
	return r
}

// MotifCount counts all occurrences of a motif in a sequence,
// including overlaps.
func MotifCount(m, s []byte) (x int) {
	for searched := 0; ; searched++ {
		i := bytes.Index(s[searched:], m)
		if i < 0 {
			break
		}
		x++
		searched += i
	}
	return
}

// ModalKmer returns the most frequent k-mers in a string.
//
//	k is the k-mer length
//	s is the string
func ModalKmer(k int, s string) (m []string) {
	f := map[string]int{}
	max := 0
	for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
		kmer := s[i:j]
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
	return
}

func KmerClumps(k, L, t int, s string) []string {
	cs := map[string]bool{}  // clump set. found clumps.
	pm := map[string][]int{} // position map. start positions by kmer.
	w := L - k               // window for start positions
	for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
		kmer := s[i:j]
		sp := append(pm[kmer], i)
		pm[kmer] = sp
		if len(sp) >= t && i-sp[len(sp)-t] <= w {
			cs[kmer] = true
		}
	}
	c := make([]string, len(cs))
	i := 0
	for kmer := range cs {
		c[i] = kmer
		i++
	}
	return c
}

// MotifSeqDist returns the minimum hamming distance from motif m to any
// same length kmer in string s.
func MotifSeqDist(m, s []byte) int {
	min := len(m)
	for i, j := 0, len(m); j < len(s); i, j = i+1, j+1 {
		if h := Hamming(m, s[i:j]); h < min {
			min = h
		}
	}
	return min
}

// MotifSetDist is a distance measure from a motif m to a set of strings l.
//
// (Not a mathematical set, just a list.)
//
// It is the sum of distances MotifSeqDist from m to each string in l.
func MotifSetDist(m []byte, l [][]byte) (d int) {
	for _, s := range l {
		d += MotifSeqDist(m, s)
	}
	return
}

// KmersNearestMotif returns the kmers in s having minimum hamming distance
// from motif m.
func KmersNearestMotif(m, s []byte) (k [][]byte) {
	min := len(m)
	for i, j := 0, len(m); j < len(s); i, j = i+1, j+1 {
		switch h := Hamming(m, s[i:j]); {
		case h < min:
			min = h
			k = [][]byte{s[i:j]}
		case h == min:
			k = append(k, s[i:j])
		}
	}
	return
}

// HammingAllIndex returns indexes of all kmers in s that are within
// hamming distance d of motif m.
func HammingAllIndex(d int, m, s []byte) (x []int) {
	for i, j := 0, len(m); j <= len(s); i, j = i+1, j+1 {
		if h := Hamming(m, s[i:j]); h <= d {
			x = append(x, i)
		}
	}
	return
}

// HammingCommonMotifs finds a list of motifs of length k that appear within
// hamming distance d in each string in the list dna.
func HammingCommonMotifs(dna []string, k, d int) (r []string) {
	// collect unique kmers
	uk := map[string]struct{}{}
	for _, s := range dna {
		for i, j := 0, k; j < len(s); i, j = i+1, j+1 {
			uk[s[i:j]] = struct{}{}
		}
	}
	// collect unique variants
	kb := make(DNA8, k)
	v := map[string]struct{}{}
	for k1 := range uk {
		copy(kb, k1)
		for _, ap := range kb.HammingVariants(d) {
			v[ap] = struct{}{}
		}
	}
	// test each variant
	a := make([]bool, len(dna))
l1:
	for ap := range v {
		for i := range a {
			a[i] = false
		}
		found := 0
		copy(kb, ap)
		for _, vv := range kb.HammingVariants(d) {
			for i, s := range dna {
				if a[i] {
					continue
				}
				if x := strings.Index(s, vv); x >= 0 {
					a[i] = true
					found++
					if found == len(dna) {
						r = append(r, ap)
						continue l1
					}
				}
			}
		}
	}
	return
}

func KmerComposition(k int, s string) []string {
	c := make([]string, len(s)-k+1)
	for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
		c[i] = s[i:j]
	}
	sort.Strings(c)
	return c
}

// ugh.  don't like any of this stuff
/*
type Adjacency map[string]map[string]struct{}

func Overlap(c []string) Adjacency {
	a := Adjacency{}
	for _, s := range c {
		suff := s[1:]
		for _, d := range c {
			if d[:len(d)-1] == suff {
				m, ok := a[s]
				if !ok {
					m = map[string]struct{}{}
				}
				m[d] = struct{}{}
				a[s] = m
			}
		}
	}
	return a
}

func DeBruijn(k int, s string) Adjacency {
	a := Adjacency{}
	for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
		pre := s[i : j-1]
		suf := s[i+1 : j]
		m, ok := a[pre]
		if !ok {
			m = map[string]struct{}{}
		}
		m[suf] = struct{}{}
		a[pre] = m
	}
	return a
}

func DeBruijnK(l []string) Adjacency {
	a := Adjacency{}
	for _, kmer := range l {
		pre := kmer[:len(kmer)-1]
		suf := kmer[1:]
		m, ok := a[pre]
		if !ok {
			m = map[string]struct{}{}
		}
		m[suf] = struct{}{}
		a[pre] = m
	}
	return a
}
*/

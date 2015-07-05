// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"bytes"
	"regexp"
	"sort"
	"strings"
)

// Seq is a general purpose byte sequence type.
//
// Bytes in a Seq are generally interpreted as separate symbols.  That is,
// Seqs do not generally hold multibyte symbols such as UTF-8 runes.
// Methods are generally case sensitive, in contrast to case insensitive types
// such as DNA8.
type Seq []byte

func (s Seq) String() string {
	return string(s)
}

const LCBit = 0x20 // Bit mask for ASCII lower case

// ToLower returns a new sequence with ASCII upper case converted to
// lower case.
//
// All other bytes are copied as-is.
func (s Seq) ToLower() Seq {
	r := make(Seq, len(s))
	for i, b := range s {
		if b >= 'A' && b <= 'Z' {
			b |= LCBit
		}
		r[i] = b
	}
	return r
}

// ToUpper returns a new sequence with ASCII lower case converted to
// upper case.
//
// All other bytes are copied as-is.
func (s Seq) ToUpper() Seq {
	r := make(Seq, len(s))
	for i, b := range s {
		if b >= 'a' && b <= 'z' {
			b &^= LCBit
		}
		r[i] = b
	}
	return r
}

// Freq counts the frequency of all symbols appearing in a sequence.
//
// The returned array contains a count of each symbol appearing in the sequence.
func (s Seq) Freq() *[256]int {
	var a [256]int
	for _, b := range s {
		a[b]++
	}
	return &a
}

// FreqMap returns counts of all symbols appearing in a sequence.
//
// The returned map contains a key for each symbol appearing in the sequence.
// If a symbol does not appear in the sequence, there will be no map key for it.
// Upper and lower case are treated as distinct.  If a sequence contains both
// 'A' and 'a' the result will contain separate counts for the two symbols.
func (s Seq) FreqMap() map[byte]int {
	m := map[byte]int{}
	for _, b := range s {
		m[b]++
	}
	return m
}

// case insensitive.  used by methods of both DNA8 and RNA8.
func baseFreq8(s []byte) (a, c, tu, g int) {
	var n [7]int
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
func (s Seq) Hamming(t Seq) int {
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
func (s Seq) AllIndex(m Seq) (x []int) {
	for searched := 0; ; {
		i := bytes.Index(s[searched:], m)
		if i < 0 {
			break
		}
		i += searched
		x = append(x, i)
		searched = i + 1
	}
	return
}

// AllCount counts all occurrences of a motif in a sequence,
// including overlaps.
//
// It is equivalent to len(s.AllIndex(m)).
func (s Seq) AllCount(m Seq) (x int) {
	for searched := 0; ; {
		i := bytes.Index(s[searched:], m)
		if i < 0 {
			break
		}
		x++
		searched += i + 1
	}
	return
}

// ModalKmers returns the most frequent k-mers in a string.
//
//	s is the string to search.
//	k is the k-mer length.
func ModalKmers(s string, k int) (m []string) {
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
	cs := map[string]struct{}{} // clump set. found clumps.
	pm := map[string][]int{}    // position map. start positions by kmer.
	w := L - k                  // window for start positions
	for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
		kmer := s[i:j]
		sp := append(pm[kmer], i)
		pm[kmer] = sp
		if len(sp) >= t && i-sp[len(sp)-t] <= w {
			cs[kmer] = struct{}{}
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

// KmersNearestMotif returns the kmers in s having minimum hamming distance
// from motif m.
func (s Seq) KmersNearestMotif(m Seq) (k []Seq) {
	min := len(m)
	for i, j := 0, len(m); j <= len(s); i, j = i+1, j+1 {
		switch h := m.Hamming(s[i:j]); {
		case h < min:
			min = h
			k = []Seq{s[i:j]}
		case h == min:
			k = append(k, s[i:j])
		}
	}
	return
}

// HammingAllIndex returns indexes of all kmers in s that are within
// hamming distance d of motif m.
func (s Seq) HammingAllIndex(d int, m Seq) (x []int) {
	for i, j := 0, len(m); j <= len(s); i, j = i+1, j+1 {
		if h := m.Hamming(s[i:j]); h <= d {
			x = append(x, i)
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

// Splice removes specified sequences (introns) from a longer sequence.
//
// Argument intr is the list of sequences to remove.
// The result is seq with intr sequences removed.
//
// The algorithm performs a simple one-time removal. It does not deal with
//  overlapping matches or matches on the spliced sequences.
func (seq Seq) Splice(intr []string) Seq {
	if len(intr) == 0 {
		return seq
	}
	q := make([]string, len(intr))
	for i, t := range intr {
		q[i] = regexp.QuoteMeta(t)
	}
	return regexp.MustCompile(strings.Join(q, "|")).ReplaceAllLiteral(seq, nil)
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

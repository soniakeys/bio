// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"bytes"
	"regexp"
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

// KmerComposition returns the frequency of kmers occurring in s.
//
// Strings are indexed by byte, not rune.
//
// Compare to DNA8.FreqArray which also computes kmer frequencies.  Most
// significantly, KmerComposition can handle large values of k.
func KmerComposition(k int, s string) map[string]int {
	c := map[string]int{}
	for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
		c[s[i:j]]++
	}
	return c
}

// KCompositionDist computes a distance metric between two strings, commonly
// called "k-tuple distance."
//
// The metric does not use an alignment of the strings but instead compares
// their kmer composition, as returned by KmerComposition.
//
// The result is the sum of absolute values of difference in frequency of
// kmers occurring in the two strings.
func KCompositionDist(k int, l, m string) (d int) {
	c := KmerComposition(k, l)
	for i, j := 0, k; j <= len(m); i, j = i+1, j+1 {
		c[m[i:j]]--
	}
	for _, n := range c {
		if n > 0 {
			d += n
		} else {
			d -= n
		}
	}
	return
}

// KCompositionDistMat computes a distance matrix for a list of strings.
//
// The distance metric is KCompositionDist (k-tuple distance).
//
// While KCompositionDist returns an integer distance, the distance matrix
// returned here uses float64s, for easy conversion to the DistanceMatrix
// type of package `cluster`.
func KCompositionDistMat(k int, l []string) [][]float64 {
	c := make([]map[string]int, len(l))
	for i, s := range l {
		c[i] = KmerComposition(k, s)
	}
	d := make([][]float64, len(l))
	d[0] = make([]float64, len(l))
	for i := 1; i < len(c); i++ {
		di := make([]float64, len(l))
		d[i] = di
		ci := c[i]
		for j, cj := range c[:i] {
			sum := 0
			for kmer, n := range ci {
				n -= cj[kmer]
				if n > 0 {
					sum += n
				} else {
					sum -= n
				}
			}
			for kmer, n := range cj {
				if _, ok := ci[kmer]; !ok {
					sum += n
				}
			}
			f := float64(sum)
			di[j] = f
			d[j][i] = f
		}
	}
	return d
}

// Splice removes specified sequences (introns) from a longer sequence.
//
// Argument intr is the list of sequences to remove.
// The result is seq with intr sequences removed.
//
// The algorithm performs a simple one-time removal. It does not deal with
// overlapping matches or matches on the spliced sequences.
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

// KMP type represents a prepocessed pattern as used by the
// Knuth-Morris-Pratt string search algorithm.
//
// Member P is the Seq passed to the constructor.
// KMP member functions will not work if the Seq is subsequently modified.
// Copy the Seq if needed.
type KMP struct {
	P Seq
	b []int
}

// NewKMP creates and generates a KMP object from pattern p.
//
// Argument P is stored as a struct member.
// KMP member functions will not work if the Seq is subsequently modified.
// Copy the Seq if needed.
//
// O(len(P)) time complexity.
func NewKMP(P Seq) KMP {
	// ref: http://www.inf.fh-flensburg.de/lang/algorithmen/pattern/kmpen.htm
	b := make([]int, len(P)+1)
	i, j := 0, -1
	b[0] = j
	for i < len(P) {
		for j >= 0 && P[i] != P[j] {
			j = b[j]
		}
		i++
		j++
		b[i] = j
	}
	return KMP{P, b}
}

// Index returns the slice index of the first match of k.P in t.
//
// It returns -1 if there is no match.
//
// O(len(t)) time complexity.
func (k KMP) Index(t Seq) int {
	i, j := 0, -1
	for i < len(t) {
		for j >= 0 && t[i] != k.P[j] {
			j = k.b[j]
		}
		i++
		j++
		if j == len(k.P) {
			return i - j
		}
	}
	return -1
}

// AllIndex returns indexes of all matches of k.P in t.
//
// O(len(t)) time complexity.
func (k KMP) AllIndex(t Seq) (x []int) {
	i, j := 0, -1
	for i < len(t) {
		for j >= 0 && t[i] != k.P[j] {
			j = k.b[j]
		}
		i++
		j++
		if j == len(k.P) {
			x = append(x, i-j)
			j = k.b[j]
		}
	}
	return
}

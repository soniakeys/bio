// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"errors"

	"github.com/soniakeys/graph"
)

// Str is a sequence type based on the string type.
//
// As with the Seq type, bytes are interpreted as separate symbols.  That is,
// Strs do not hold multibyte symbols such as UTF-8 runes.  Also as with the
// Seq type, methods are generally case sensitive.
//
// Often a method is written to the str type rather than a byte-based sequence
// type because the implmentation uses maps.
type Str string

// StrList is a list of Strs, not necessarily of the same length.
//
// See StrKmers for Strs of all the same length.
type StrList []Str

// StrKmers is a list of Strs where each is expected to be the same length (k).
//
// See StrList for Strs that may have varied lengths.
type StrKmers []Str

// StrFreq is a frequency map of Strs.
//
// The usual meaning is that for map m and str s, m[s] is the occurrence
// frequency of s in some context.
type StrFreq map[Str]int

// Hamming returns the Hamming distance between two Strs.
// Comparison is done byte-wise and so is case sensitive.
//
// Strs s and t must be of equal length.
// Panic or nonsense results if the lengths are unequal.
func (s Str) Hamming(t Str) int {
	h := 0
	for i, b := range []byte(s) {
		if b != t[i] {
			h++
		}
	}
	return h
}

// DNA8HammingVariants computes DNA k-mers -- as strings -- within hamming
// distance d of receiver kmer k.
//
// Argument kmer must have DNA8 content.
// Variants are appended to argument v and returned.
// Like the the DNA8 version of the function, case is preserved by position.
// Note though that depending on your application, you may need to supply a
// kmer of all the same case so that you get variants of all the same case.
//
// See also DNA8.HammingVariants.
func (kmer Str) DNA8HammingVariants(d int, v StrKmers) StrKmers {
	// recursive, but minimizes allocations
	v = append(v, kmer)
	if d == 0 {
		return v
	}
	k8 := DNA8(kmer)
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
				v = append(v, Str(k8))
				if h > 1 && len(sub) > 1 {
					f(sub[1:], h-1)
				}
				vb += 2
			}
			sub[0] = b
		}
	}
	f(k8, d)
	return v
}

// ModalKmers returns the most frequent k-mers in a string.
//
// s is the string to search, k is the k-mer length.
func (s Str) ModalKmers(k int) (m StrKmers) {
	f := StrFreq{}
	max := 0
	for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
		kmer := s[i:j]
		n := f[kmer] + 1
		f[kmer] = n
		switch {
		case n > max:
			m = StrKmers{kmer}
			max = n
		case n == max:
			m = append(m, kmer)
		}
	}
	return
}

// ProximalKmerRepeats finds kmers that repeat within some proximity of each
// other.
//
// Argument k is the kmer length, L is the window of proximity in which repeats
// must be found, and t is the number of repeats that must be found within a
// window.
//
// Returned is a list of kmers found to repeat at least t time within some
// window of length L somewhere in receiver string s.  (Sorry, no position
// information is returned currently.)
func (s Str) ProximalKmerRepeats(k, L, t int) StrKmers {
	cs := map[Str]struct{}{} // clump set. found clumps.
	pm := map[Str][]int{}    // position map. start positions by kmer.
	w := L - k               // window for start positions
	for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
		kmer := s[i:j]
		sp := append(pm[kmer], i)
		pm[kmer] = sp
		if len(sp) >= t && i-sp[len(sp)-t] <= w {
			cs[kmer] = struct{}{}
		}
	}
	c := make(StrKmers, len(cs))
	i := 0
	for kmer := range cs {
		c[i] = kmer
		i++
	}
	return c
}

// KmerComposition returns the frequency of kmers occurring in s.
//
// Strings are indexed by byte, not rune.
//
// Compare to DNA8.FreqArray which also computes kmer frequencies.  Most
// significantly, KmerComposition can handle large values of k.
func (s Str) KmerComposition(k int) StrFreq {
	c := StrFreq{}
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
func (l Str) KCompositionDist(k int, m Str) (d int) {
	c := l.KmerComposition(k)
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
func (l StrList) KCompositionDistMat(k int) [][]float64 {
	c := make([]StrFreq, len(l))
	for i, s := range l {
		c[i] = s.KmerComposition(k)
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

// OverlapKmers constructs a sequence from overlapping kmers.
//
// Argument `order` contains indexes into `kmers` and determines a
// sequence in which kmers overlap, offset by 1.
//
// The function does not validate overlap consistency.
func (kmers StrKmers) OverlapKmers(order []int) (Seq, error) {
	if len(order) == 0 {
		return nil, nil
	}
	if len(kmers) == 0 {
		return nil, errors.New("empty kmer list")
	}
	k := len(kmers[0])
	if k == 0 {
		return nil, errors.New("zero length kmers")
	}
	for _, m := range kmers[1:] {
		if len(m) != k {
			return nil, errors.New("kmers not all the same length")
		}
	}
	last := len(order) - 1
	s := make(Seq, last+k)
	for i, o := range order[:last] {
		s[i] = kmers[o][0]
	}
	copy(s[last:], kmers[order[last]])
	return s, nil
}

type deBruijn struct {
	jmers StrKmers
	jNode map[Str]int
	graph [][]int
}

func newDeBruijn() *deBruijn {
	return &deBruijn{jNode: map[Str]int{}}
}

func (d *deBruijn) node(jmer Str) int {
	if n, ok := d.jNode[jmer]; ok {
		return n
	}
	n := len(d.jmers)
	d.jNode[jmer] = n
	d.jmers = append(d.jmers, jmer)
	d.graph = append(d.graph, nil)
	return n
}

func (d *deBruijn) arc(fr, to int) {
	d.graph[fr] = append(d.graph[fr], to)
}

// DeBruijn constructs a DeBruijn graph from a multiset
// of kmers, for example the kmer composition of a string.
//
// The graph represents overlaps of length k-1 of the kmers.
// Nodes of the returned graph are labeled by
// k-1-mer (jmer) substrings of the kmers.
func (kmerFreq StrFreq) DeBruijn() (graph [][]int, jmers StrKmers, err error) {
	if len(kmerFreq) == 0 {
		return
	}
	var k int
	for kmer := range kmerFreq {
		k = len(kmer)
		break
	}
	j := k - 1
	d := newDeBruijn()
	for kmer, mult := range kmerFreq {
		if len(kmer) != k {
			return nil, nil, errors.New("kmers have different length")
		}
		fr := d.node(kmer[:j])
		to := d.node(kmer[1:])
		for n := 0; n < mult; n++ {
			d.arc(fr, to)
		}
	}
	return d.graph, d.jmers, nil
}

// DeBruijn constructs a DeBruijn graph from a list of kmers.
//
// The graph represents overlaps of length k-1 of the kmers.
// Nodes of the returned graph are labeled by
// k-1-mer (jmer) substrings of the kmers.
func (kmers StrKmers) DeBruijn() (graph [][]int, jmers StrKmers, err error) {
	if len(kmers) == 0 {
		return
	}
	k := len(kmers[0])
	j := k - 1
	d := newDeBruijn()
	for _, kmer := range kmers {
		if len(kmer) != k {
			return nil, nil, errors.New("kmers have different length")
		}
		fr := d.node(kmer[:j])
		to := d.node(kmer[1:])
		d.arc(fr, to)
	}
	return d.graph, d.jmers, nil
}

// DeBruijn constructs a DeBruijn graph from the kmer composition
// of receiver string s.
//
// The graph represents overlaps of the sequence of kmers in the string.
// Overlaps are of length k-1.  Nodes of the returned graph are labeled by
// k-1-mer (jmer) substrings of the kmers.
func (s Str) DeBruijn(k int) (graph [][]int, jmers StrKmers) {
	d := newDeBruijn()
	to := d.node(s[:k-1])
	for i, j := 1, k; j <= len(s); i, j = i+1, j+1 {
		fr := to
		to = d.node(s[i:j])
		d.arc(fr, to)
	}
	return d.graph, d.jmers
}

// Contigs finds contigs in the DeBruijn graph represented by kmer list kmers.
func (kmers StrKmers) Contigs() ([]Seq, error) {
	var g graph.AdjacencyList
	g, jmers, err := kmers.DeBruijn()
	if err != nil {
		return nil, err
	}
	ps := g.MaximalNonBranchingPaths()
	s := make([]Seq, len(ps))
	for i, p := range ps {
		s[i], _ = jmers.OverlapKmers(p) // (length already checked)
	}
	return s, nil
}

type DistFuncStr func(Str, Str) float64

// DistanceMatrix computes a distance matrix for a StrList and a DistFuncStr that
// compares two Strs.
func (l StrList) DistanceMatrix(f DistFuncStr) [][]float64 {
	// code identical to DNA8List.DistanceMatrix
	d := make([][]float64, len(l))
	d[0] = make([]float64, len(l))
	for i := 1; i < len(l); i++ {
		li := l[i]
		di := make([]float64, len(l))
		d[i] = di
		for j, lj := range l[:i] {
			dij := f(li, lj)
			di[j] = dij
			d[j][i] = dij
		}
	}
	return d
}

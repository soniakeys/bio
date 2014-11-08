// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"bytes"
	"math"
	"math/rand"
	"sort"
)

// DNA8List is a general purpose list of DNA8 strings.
type DNA8List []DNA8

// Kmers is a list of DNA8 strings all the same length.
//
// Methods on the type may require a non-empty list.  Panic or nonsense may
// result from calling methods on an empty list or a list with unequal length
// strings.  Test with method Uniform where this is important.
type Kmers []DNA8

// Uniform returns true if receiver list x is is non-empty and contains strings
// all the same length.
func (x Kmers) Uniform() bool {
	if len(x) == 0 {
		return false
	}
	s0 := x[0]
	for _, s := range x[1:] {
		if len(s) != len(s0) {
			return false
		}
	}
	return true
}

// ConsensusHamming is the cumulative sum of hamming distances between the
// kmers in x and a consensus string formed from them.
//
// Alternatively described, it is the number of unconserved bases.
// A base is conserved when it is the modal base in a position.
// If multimodal, only a single base is considered conserved, others are
// counted as unconserved.
//
// It represents a measure how much difference exists among matrix of uniform
// length strings.
func (x Kmers) ConsensusHamming() (s int) {
	// algorithm by alternative description, "Score"
	k := len(x[0])
	for i := 0; i < k; i++ {
		var f [7]int
		max := 0
		for _, m := range x {
			b2 := m[i] & 6
			f[b2]++
			if f[b2] > max {
				max = f[b2]
			}
		}
		s += len(x) - max
	}
	return
}

// CountProfile represents base counts by position relative some set of
// DNA strings.
//
// The nominal base order is ACTG.
type CountProfile [][4]int

// CountProfile builds a CountProfile object from Kmers x.
func (x Kmers) CountProfile() CountProfile {
	// compare to DNA8List.Profile()
	c := make(CountProfile, len(x[0]))
	for _, m := range x {
		for i, b := range m {
			c[i][b>>1&3]++
		}
	}
	return c
}

// Add adds DNA sequence s to the receiver CountProfile p.  Mixed case is
// allowed.  String lengths not equal to the profile matrix length are allowed,
// with shortages and excesses ignored.
// Symbols other than DNA bases are ignored.
func (p CountProfile) Add(s DNA) {
	if len(s) > len(p) {
		s = s[:len(p)]
	}
	for x, b := range s {
		switch lb := b | LCBit; lb {
		case 'a', 'c', 't', 'g':
			p[x][lb>>1&3]++
		}
	}
}

// Consensus returns a consensus string and difference score from a populated
// CountProfile.
// If all base counts are zero in a particular position,
// a '-' is emitted in that position.
func (pm CountProfile) Consensus() (cs DNA, score int) {
	cs = make(DNA, len(pm))
	for x := range cs {
		c := byte('-') // consensus symbol
		fc := 0        // freq of consensus symbol
		if f := pm[x][0]; f > fc {
			c = 'A'
			fc = f
		}
		if f := pm[x][1]; f > fc {
			c = 'C'
			fc = f
		}
		if f := pm[x][2]; f > fc {
			c = 'T'
			fc = f
		}
		if f := pm[x][3]; f > fc {
			c = 'G'
			fc = f
		}
		cs[x] = c
		score += fc
	}
	return
}

// FracProfile represents a profile matrix corresponding to a list of DNA
// k-mers.
//
// Elements values correspond to the fraction of each base by position.
// Row order is ACTG.
type FracProfile [][4]float64

// FracProfile constructs the FracProfile matrix corresponding to receiver
// Kmers x.
//
// In the resulting profile, element values for the four bases will sum to
// approximately 1.0 for each position.
func (x Kmers) FracProfile() FracProfile {
	// compare to DNA8List.BaseCount
	p := make(FracProfile, len(x[0]))
	inc := 1 / float64(len(x))
	for _, m := range x {
		for i, b := range m {
			p[i][b>>1&3] += inc
		}
	}
	return p
}

// LaplaceProfile constructs the profile matrix corresponding to
// motif matrix m, augmented by Laplace's Rule of Succession.  That is,
// with a pseudocount of 1 added to each base count.
func (x Kmers) LaplaceProfile() FracProfile {
	inc := 1 / float64(len(x)+4)
	p := NewPseudoProfile(len(x[0]), inc)
	for _, m := range x {
		for i, b := range m {
			p[i][b>>1&3] += inc
		}
	}
	return p
}

// NewPseudoProfile allocates a FracProfile matrix for string length k,
// initializing all elements with pseudocount probabilty p.
func NewPseudoProfile(k int, p float64) FracProfile {
	a := make(FracProfile, k)
	p4 := [4]float64{p, p, p, p}
	for i := range a {
		a[i] = p4
	}
	return a
}

// KmerProb computes the probability of a kmer given profile p.
func (p FracProfile) KmerProb(kmer DNA8) float64 {
	pr := 1.
	for i, b := range kmer {
		pr *= p[i][b>>1&3]
	}
	return pr
}

// ProbKmer returns the profile-most probable kmer in s.
func (p FracProfile) ProbKmer(s DNA8) (k DNA8) {
	max := 0.
	i := 0
	j := len(p)
	for j <= len(s) {
		if pr := p.KmerProb(s[i:j]); pr > max {
			k = s[i:j]
			max = pr
		}
		i++
		j++
	}
	return
}

// RandWeighted returns an int >= 0 and < len(weights)
//
// The probability of returning a value n is proportional to weight[n].
func RandWeighted(weights []float64) (n int) {
	sum := 0.
	for _, w := range weights {
		sum += w
	}
	f := rand.Float64() * sum
	c := 0.
	for v, w := range weights {
		c += w
		if c > f {
			return v
		}
	}
	return len(weights) - 1
}

func (p FracProfile) RandKmer(s DNA8) DNA8 {
	k := len(p[0])
	c := make([]float64, len(s)-k+1)
	c[0] = p.KmerProb(s[:k])
	for i, j := 1, k+1; j < len(s); i, j = i+1, j+1 {
		c[i] = c[i-1] + p.KmerProb(s[i:j])
	}
	i := sort.SearchFloat64s(c, c[len(c)-1]*rand.Float64())
	return s[i : i+k]
}

func GibbsSampler(l []DNA8, k, N int) (Kmers, int) {
	motifs := make(Kmers, len(l))
	for i, s := range l {
		j := rand.Intn(len(s) - k + 1)
		motifs[i] = s[j : j+k]
	}
	best := len(l) * k
	bestMotifs := make(Kmers, len(l))
	for i := 0; i < N; i++ {
		j := rand.Intn(len(l))
		motifs[j] = motifs[0]
		motifs[0], motifs[j] = motifs[j], motifs[1:].FracProfile().RandKmer(l[j])
		if s := motifs.ConsensusHamming(); s < best {
			best = s
			copy(bestMotifs, motifs)
		}
	}
	return bestMotifs, best
}

func (p FracProfile) Motifs(l []DNA8) Kmers {
	m := make(Kmers, len(l))
	for i, s := range l {
		m[i] = p.ProbKmer(s)
	}
	return m
}

func (l DNA8List) RandomMotifSearch(k int) (Kmers, int) {
	motifs := make(Kmers, len(l))
	for i, s := range l {
		j := rand.Intn(len(s) - k + 1)
		motifs[i] = s[j : j+k]
	}
	bestMotifs := make(Kmers, len(l))
	best := len(l) * k
	for {
		motifs = motifs.LaplaceProfile().Motifs(l)
		if s := motifs.ConsensusHamming(); s < best {
			copy(bestMotifs, motifs)
			best = s
		} else {
			return bestMotifs, s
		}
	}
}

// Consensus generates a consensus string from receiver Kmers c.
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
func (c Kmers) Consensus() (seq DNA8, score int) {
	s := c[0]
	r := make(DNA8, len(s))
	// compute position by position, without constructing profile matrix
	for i := range r {
		// profile single position
		var n [7]int
		for _, s := range c {
			if i < len(s) {
				// (compare to DNA version)
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
		// (compare to DNA version)
		r[i] = "A C T G"[maxb]
		score += max
	}
	return r, score
}

// Entropy computes entropy for a motif matrix.  It is the sum of entropies
// in each column.
func (x Kmers) Entropy() (e float64) {
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

// Median returns a list of kmers that are at minimum distance
// (by the method MotifHamming) to a list of strings l.
//
// Sequences in receiver l may be of different lengths.
//
// The algorithm is brute force and practical only when k is small.
func (l Kmers) Median(k int) (m Kmers) {
	z := DNA8(bytes.Repeat([]byte{'A'}, k))
	min := len(l[0])
	for p := append(DNA8{}, z...); ; {
		switch d := DNA8List(l).MotifHamming(p); {
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
func (l DNA8List) MotifSearch(k int) (m Kmers) {
	bestMotifs := make(Kmers, len(l))
	for i, s := range l {
		bestMotifs[i] = s[:k]
	}
	bestScore := bestMotifs.ConsensusHamming()
	motifs := make(Kmers, len(l))
	s0 := l[0]
	i := 0
	j := k
	for j <= len(s0) {
		motifs[0] = s0[i:j]
		for i := 1; i < len(l); i++ {
			p := motifs[:i].LaplaceProfile()
			motifs[i] = p.ProbKmer(l[i])
		}
		if s := motifs.ConsensusHamming(); s < bestScore {
			bestScore = s
			copy(bestMotifs, motifs)
		}
		i++
		j++
	}
	return bestMotifs
}

// DNAConsensus returns a consensus sequence from multiple sequences.
//
// Consensus in each position is the most frequent base in that
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
	// compute position by position, without constructing profile matrix
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

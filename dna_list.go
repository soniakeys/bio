// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"bytes"
	"fmt"
	"math"
	"math/rand"
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

//----------------------------------------------------------------------------
// Profile and consensus synopsis
//
// DNAConsensus function,  works on []DNA, computes consensus and "match" score
//                         directly from strings without constructing profile.
// Kmers.Consensus,        as above, but for DNA8.
// Kmers.ConsensusHamming, computes cumulative hamming distance ("mismatch")
//                         directly from kmers without constructing profile.
// Kmers.Entropy,          computes Entropy--"mismatch" sense like hamming
//                         but better.
// Kmers.EntropyContributions, individual entropy terms.
// FracProfile.Entropy
// FracProfile.CrossEntropy
// FracProfile.RelativeEntropy
//
// Kmers.CountProfile,     constructs and populates from DNA8.
// CountProfile.Add,       adds a DNA (not DNA8) string.
// CountProfile.Consensus, returns DNA consensus and "match" score.
//
// Kmers.FracProfile,      constructs and populates from DNA8.
// Kmers.LaplaceProfile,   same as FracProfile, but with pseudocounts.
// FracProfile.KmerProbability(kmer DNA8) float64
// FracProfile.MostProbKmer(s DNA8) (kmer DNA8)
// FracProfile.MostProbKmers(l []DNA8) Kmers, as above, just broadcast
//----------------------------------------------------------------------------
// Motif finding synopsis
//
// DNA8List.HammingMotifs  brute force over kmers present
// DNA8List.HammingMotifs1 brute force over kmers present
// DNA8List.HammingMotifs2 brute force over kmers present
// DNA8List.MedianMotifs,  brute force over all possible kmers
// DNA8List.GreedyMotifSearch
// DNA8List.RandomMotifSearch
// DNA8List.GibbsMotifSearch
//----------------------------------------------------------------------------

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

// ConsensusHamming is the cumulative sum of hamming distances between the
// kmers in x and a consensus string formed from them.
//
// Alternatively described, it is the number of unconserved bases.
// A base is conserved when it is the consenus, or modal, base in a position.
// If multimodal, only a single base is considered conserved, others are
// counted as unconserved.
//
// Either way, it represents a measure how much difference exists among a
// matrix of uniform length strings.
func (x Kmers) ConsensusHamming() (s int) {
	// algorithm by alternative description
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

// EntropyContributions computes contributions of each base in each position
// to the total entropy for a motif matrix.
func (x Kmers) EntropyContributions() [][4]float64 {
	k := len(x[0])
	e := make([][4]float64, k)
	inc := 1 / float64(len(x))
	for i := range e {
		var f [4]float64
		for _, m := range x {
			f[m[i]>>1&3] += inc
		}
		for j, p := range f {
			if p > 0 {
				e[i][j] = -p * math.Log2(p)
			}
		}
	}
	return e
}

// CountProfile represents base counts by position among a set of
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

// Consensus returns a consensus string and count from a populated CountProfile.
//
// If all base counts are zero in a particular position,
// a '-' is emitted in that position in the consensus result.
//
// The result 'count' is the sum of consensus counts.  (Note this is the
// opposite sense from that returned by Kmers.ConsensusHamming.  This returns
// a measure of match, Kmers.ConsensusHamming returns a measure of mismatch.)
func (pm CountProfile) Consensus() (cs DNA, count int) {
	cs = make(DNA, len(pm))
	for x := range cs {
		c := byte('-') // consensus symbol
		fc := 0        // freq of consensus symbol
		for bx, f := range pm[x] {
			if f > fc {
				c = "ACTG"[bx]
				fc = f
			}
		}
		cs[x] = c
		count += fc
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
	p := newPseudoProfile(len(x[0]), inc)
	for _, m := range x {
		for i, b := range m {
			p[i][b>>1&3] += inc
		}
	}
	return p
}

// newPseudoProfile allocates a FracProfile matrix for string length k,
// initializing all elements with pseudocount probabilty p.
func newPseudoProfile(k int, p float64) FracProfile {
	a := make(FracProfile, k)
	p4 := [4]float64{p, p, p, p}
	for i := range a {
		a[i] = p4
	}
	return a
}

// Entropy computes entropy for a FracProfile.  It is the sum of entropies
// in each column.
func (p FracProfile) Entropy() (e float64) {
	for i := range p {
		for _, pr := range p[i] {
			if pr > 0 {
				e -= pr * math.Log2(pr)
			}
		}
	}
	return
}

// BaseFreq returns fractional base frequencies in receiver list l.
//
// The four elements of the result sum to 1.0.
func (l DNA8List) BaseFreq() (b [4]float64) {
	var n [7]int
	for _, s := range l {
		for _, b := range s {
			n[b&6]++
		}
	}
	c := 1 / float64(n[0]+n[2]+n[4]+n[6])
	for i := range b {
		b[i] = float64(n[i*2]) * c
	}
	return
}

// CrossEntropy computes cross-entropy for a FracProfile.
func (p FracProfile) CrossEntropy(log2b [4]float64) (e float64) {
	for i := range p {
		for j, pr := range p[i] {
			if pr > 0 {
				e -= pr * log2b[j]
			}
		}
	}
	return
}

// RelativeEntropy computes cross-entropy for a FracProfile.
func (p FracProfile) RelativeEntropy(b [4]float64) (e float64) {
	for i := range p {
		for j, pr := range p[i] {
			if pr > 0 {
				e += pr * math.Log2(pr/b[j])
			}
		}
	}
	return
}

// KmerProbability computes the probability of a kmer given profile p.
func (p FracProfile) KmerProbability(kmer DNA8) float64 {
	pr := 1.
	for i, b := range kmer {
		pr *= p[i][b>>1&3]
	}
	return pr
}

// MostProbKmer returns the profile-most probable kmer in s.
//
// If string s is shorter than profile p, MostProbKmer returns nil.
func (p FracProfile) MostProbKmer(s DNA8) (kmer DNA8) {
	max := -1.
	end := len(s) - len(p)
	for i := 0; i <= end; i++ {
		k := s[i : i+len(p)]
		if pr := p.KmerProbability(k); pr > max {
			kmer = k
			max = pr
		}
	}
	return
}

func (p FracProfile) MostProbKmers(l []DNA8) Kmers {
	m := make(Kmers, len(l))
	for i, s := range l {
		m[i] = p.MostProbKmer(s)
	}
	return m
}

// HammingMotifs finds a list of motifs of length k that appear within
// hamming distance d in each string in samples.
//
// That is, for each result kmer in r, there will be a kmer within
// hamming distance d in each string of samples.
//
// Reference Compeau 2014, p. 87, Algorithm "MotifEnumeration".
func (samples DNA8List) HammingMotifs(k, d int) []DNA8 {
	// by the book...
	Patterns := map[string]DNA8{}
	for _, seq := range samples {
		for i, j := 0, k; j <= len(seq); i, j = i+1, j+1 {
		hv:
			for _, Patternʹ := range seq[i:j].HammingVariantsRef(d) {
				for _, seq := range samples {
					if h := seq.MotifHamming(Patternʹ); h > d {
						continue hv
					}
				}
				Patterns[string(Patternʹ)] = Patternʹ
			}
		}
	}
	r := make([]DNA8, len(Patterns))
	i := 0
	for _, r[i] = range Patterns {
		i++
	}
	return r
}

// HammingMotifs1.  Same result, different algorithm.
func (samples DNA8List) HammingMotifs1(k, d int) (r []string) {
	// collect unique kmers
	uk := map[string]struct{}{}
	for _, d := range samples {
		s := string(d)
		for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
			uk[s[i:j]] = struct{}{}
		}
	}
	// collect unique variants
	kb := make(DNA8, k)
	v := map[string]struct{}{}
	for k1 := range uk {
		copy(kb, k1)
		for _, ap := range kb.HammingVariants(d) {
			v[string(ap)] = struct{}{}
		}
	}
	// test each variant
	a := make([]bool, len(samples))
l1:
	for ap := range v {
		for i := range a {
			a[i] = false
		}
		found := 0
		copy(kb, ap)
		for _, vv := range kb.HammingVariants(d) {
			for i, s := range samples {
				if a[i] {
					continue
				}
				if x := bytes.Index(s, vv); x >= 0 {
					a[i] = true
					found++
					if found == len(samples) {
						r = append(r, ap)
						continue l1
					}
				}
			}
		}
	}
	return
}

// HammingMotifs2.  Same result, different algorithm.
func (samples DNA8List) HammingMotifs2(k, d int) (r []string) {
	type kset map[string]struct{}
	// dna neigbors, list parallel to dna
	// dnbs[i] is set of unique neighbors over all kmers of dna[i]
	dnbs := make([]kset, len(samples))
	// unbs, set of all unique neighbors of kmers of dna.
	unbs := kset{}

	// build dnbs, unbs in one pass.
	//
	// useful in the process is unique kmer neighbors:
	// key is a kmer of a seq of dna
	// val is neighbors of kmer
	uknbs := map[string][]string{}
	for n, seq := range samples {
		fmt.Println("seq", n, ":", seq[:10], "...")
		// unique kmers over the current sequence
		uk := kset{}
		// unique neigbors over the current sequence
		ds := kset{}
		// kmers in seq
		s := string(seq)
		for i, j := 0, k; j <= len(s); i, j = i+1, j+1 {
			kmer := s[i:j]
			if _, ok := uk[kmer]; ok {
				continue
			}
			fmt.Println("   new kmer this seq:", kmer)
			// first time kmer seen this seq
			uk[kmer] = struct{}{}
			nbs, ok := uknbs[kmer]
			if !ok {
				// first time kmer seen in dna
				vs := seq[i:j].HammingVariants(d)
				nbs = make([]string, len(vs))
				for i, v := range vs {
					nb := string(v)
					unbs[nb] = struct{}{} // store in neighbors of dna
					nbs[i] = nb           // accumulate kmer neigbors
				}
				uknbs[kmer] = nbs // store neighbors for kmer
			}
			// now whether it was new overall or not, it was
			// new to this seq so store in neighbors of this seq.
			for _, nb := range nbs {
				ds[nb] = struct{}{}
			}
		}
		// all neighbors of seq now known.
		dnbs[n] = ds
	}
	// dnbs, unbs complete.
	// test each unique neighbor
u:
	for u := range unbs {
		for _, nbs := range dnbs {
			if _, ok := nbs[u]; !ok {
				continue u
			}
		}
		r = append(r, u)
	}
	return
}

// MedianMotifs returns a list of kmers that are at minimum cumulative distance
// (by the method DNA8List.MotifHamming) to a list of strings l.
//
// Returned is a list of all motifs found with minimum MotifHamming.
// Also returned is the corresponding minimum MotifHamming cumulative distance.
//
// The algorithm is brute force and practical only when k is small.
func (l DNA8List) MedianMotifs(k int) (m Kmers, hamming int) {
	hamming = k * len(l)
	p0 := l[0][:k]
	for p := append(DNA8{}, p0...); ; {
		switch h := l.MotifHamming(p); {
		case h < hamming:
			hamming = h
			m = Kmers{append(DNA8{}, p...)}
		case h == hamming:
			m = append(m, append(DNA8{}, p...))
		default:
		}
		p.Inc()
		if bytes.Equal(p, p0) {
			break
		}
	}
	return
}

func (l DNA8List) GreedyMotifSearch(k int) (m Kmers) {
	bestMotifs := make(Kmers, len(l))
	for i, s := range l {
		bestMotifs[i] = s[:k]
	}
	bestScore := bestMotifs.ConsensusHamming()
	motifs := make(Kmers, len(l))
	s0 := l[0]
	end := len(s0) - k
	for i := 0; i <= end; i++ {
		motifs[0] = s0[i : i+k]
		for i := 1; i < len(l); i++ {
			p := motifs[:i].LaplaceProfile()
			motifs[i] = p.MostProbKmer(l[i])
		}
		if s := motifs.ConsensusHamming(); s < bestScore {
			bestScore = s
			copy(bestMotifs, motifs)
		}
	}
	return bestMotifs
}

func (l DNA8List) RandomKmers(k int) Kmers {
	kmers := make(Kmers, len(l))
	for i, s := range l {
		j := rand.Intn(len(s) - k + 1)
		kmers[i] = s[j : j+k]
	}
	return kmers
}

func (l DNA8List) ConvergedRandomMotifs(k int) (motifs Kmers, hamming int) {
	motifs = l.RandomKmers(k)
	hamming = motifs.ConsensusHamming()
	// then converge
	for {
		m := motifs.LaplaceProfile().MostProbKmers(l)
		h := m.ConsensusHamming()
		if h == hamming {
			return
		}
		motifs = m
		hamming = h
	}
}

func (l DNA8List) RandomMotifSearch(k, N int) (motifs Kmers) {
	min := k * len(l)
	for i := 0; i < N; i++ {
		m, h := l.ConvergedRandomMotifs(k)
		if h < min {
			min = h
			motifs = m
		}
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
		if c >= f {
			return v
		}
	}
	return len(weights) - 1
}

func (l DNA8List) MaxLen() (m int) {
	for _, s := range l {
		if len(s) > m {
			m = len(s)
		}
	}
	return
}

func (l DNA8List) GibbsSampler(k, N int) (motifs Kmers, hamming int) {
	motifs = l.RandomKmers(k)
	hamming = motifs.ConsensusHamming()
	m := append(Kmers{}, motifs...)
	p := make([]float64, l.MaxLen()-k+1)
	for j := 0; j < N; j++ {
		i := rand.Intn(len(l))
		m[i] = m[0]
		pf := m[1:].LaplaceProfile()
		s := l[i]
		p = p[:len(s)-k+1]
		for x := range p {
			p[x] = pf.KmerProbability(s[x : x+k])
		}
		x := RandWeighted(p)
		m[i] = s[x : x+k]
		if h := m.ConsensusHamming(); h < hamming {
			hamming = h
			copy(motifs, m)
		}
	}
	return
}

func (l DNA8List) GibbsMotifSearch(k, N, M int) (motifs Kmers) {
	min := k * len(l)
	for i := 0; i < M; i++ {
		m, h := l.GibbsSampler(k, N)
		if h < min {
			min = h
			motifs = m
		}
	}
	return
}

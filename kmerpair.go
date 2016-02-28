// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"errors"
	"fmt"

	"github.com/soniakeys/graph"
)

// KmerPair represents a kmer read pair where kmers are separated by some
// distance d.  The distance d is not represented here.  For pair p, p.A is
// the kmer closer to the beginning of some string, p.B is the kmer that
// follows after distance d.
type KmerPair struct{ A, B Str }

// ReadPairComposition generates the ideal read pair composition of a string.
//
// Argument k is the kmer size, d is the distance between paired kmers.
func (s Str) ReadPairComposition(k, d int) ReadPairFreq {
	m := map[KmerPair]int{}
	for i, jb := 0, k+d+k; jb <= len(s); i, jb = i+1, jb+1 {
		m[KmerPair{
			s[i : i+k],
			s[jb-k : jb],
		}]++
	}
	return ReadPairFreq{d, m}
}

type ReadPairList struct {
	D     int // distance between reads
	Pairs []KmerPair
}

type ReadPairFreq struct {
	D    int // distance between reads
	Freq map[KmerPair]int
}

func (pairs ReadPairList) ReadBreak(bk int) (ReadPairFreq, error) {
	ps := pairs.Pairs
	if len(ps) == 0 {
		return ReadPairFreq{}, errors.New("no read pairs")
	}
	k := len(ps[0].A)
	if bk > k {
		return ReadPairFreq{}, errors.New("bk > read k")
	}
	m := map[KmerPair]int{}
	for _, p := range ps {
		for j := 0; j <= k-bk; j++ {
			m[KmerPair{p.A[j : j+bk], p.B[j : j+bk]}]++
		}
	}
	return ReadPairFreq{pairs.D + k - bk, m}, nil
}

func (pairs ReadPairList) OverlapSeq(order []graph.NI) (Seq, error) {
	d := pairs.D
	ps := pairs.Pairs
	if len(order) <= d {
		return nil, fmt.Errorf("ordering too short (%d)", len(order))
	}
	if len(ps) == 0 {
		return nil, errors.New("empty pair list")
	}
	k := len(ps[0].A)
	fmt.Println("Overlap: ", len(order), "in order of", len(ps), "pairs.  k, d =", k, d)
	if k == 0 {
		return nil, errors.New("zero length kmers")
	}
	for _, p := range ps[1:] {
		if len(p.A) != k || len(p.B) != k {
			return nil, errors.New("kmers not all the same length")
		}
	}
	if d < 0 {
		return nil, errors.New("negative d")
	}
	last := len(order) - 1
	s := make(Seq, last+k+d+k)
	for i, o := range order[:last] {
		s[i] = ps[o].A[0]
	}
	copy(s[last:], ps[order[last]].B)
	d0 := last + k
	for i, o := range order[last-d : last] {
		s[d0+i] = ps[o].B[0]
	}
	copy(s[d0+d:], ps[order[last]].B)
	return s, nil
}

type pairDeBruijn struct {
	jpairs []KmerPair
	jNode  map[KmerPair]graph.NI
	g      graph.Directed
}

func newPairDeBruijn() *pairDeBruijn {
	return &pairDeBruijn{jNode: map[KmerPair]graph.NI{}}
}

func (pd *pairDeBruijn) node(jpair KmerPair) graph.NI {
	if n, ok := pd.jNode[jpair]; ok {
		return n
	}
	n := graph.NI(len(pd.jpairs))
	pd.jNode[jpair] = n
	pd.jpairs = append(pd.jpairs, jpair)
	pd.g.AdjacencyList = append(pd.g.AdjacencyList, nil)
	return n
}
func (pd *pairDeBruijn) arc(fr, to graph.NI) {
	pd.g.AdjacencyList[fr] = append(pd.g.AdjacencyList[fr], to)
}

// DeBruijn constructs a DeBruijn graph from a list of kmer pairs.
//
// The graph represents overlaps of length k-1 of the kmers.
// Nodes of the returned graph are labeled by
// k-1-mer (jmer) substring pairs of the kmer pairs.
func (pairs ReadPairList) DeBruijn() (g graph.Directed, jpairs ReadPairList, err error) {
	ps := pairs.Pairs
	if len(ps) == 0 {
		return
	}
	k := len(ps[0].A)
	j := k - 1
	pd := newPairDeBruijn()
	for _, pair := range ps {
		if len(pair.A) != k || len(pair.B) != k {
			err = errors.New("kmers have different length")
			return
		}
		fr := pd.node(KmerPair{pair.A[:j], pair.B[:j]})
		to := pd.node(KmerPair{pair.A[1:], pair.B[1:]})
		pd.arc(fr, to)
	}
	return pd.g, ReadPairList{pairs.D + 1, pd.jpairs}, nil
}

func (pairs ReadPairList) Contigs() (cs []Seq, err error) {
	var g graph.Directed
	g, jpairs, err := pairs.DeBruijn()
	if err != nil {
		return nil, err
	}
	var s Seq
	g.MaximalNonBranchingPaths(func(p []graph.NI) bool {
		if s, err = jpairs.OverlapSeq(p); err == nil {
			cs = append(cs, s)
			return true
		}
		return false
	})
	return
}

func (freq ReadPairFreq) Contigs() (cs []Seq, err error) {
	var g graph.Directed
	g, jpairs, err := freq.DeBruijn()
	if err != nil {
		return nil, err
	}
	var s Seq
	g.MaximalNonBranchingPaths(func(p []graph.NI) bool {
		if s, err = jpairs.OverlapSeq(p); err == nil {
			cs = append(cs, s)
			return true
		}
		return false
	})
	return
}

// DeBruijn constructs a DeBruijn graph from a frequency map of kmer pairs.
//
// The graph represents overlaps of length k-1 of the kmers.
// Nodes of the returned graph are labeled by
// k-1-mer (jmer) substring pairs of the kmer pairs.
func (freq ReadPairFreq) DeBruijn() (g graph.Directed, jpairs ReadPairList, err error) {
	m := freq.Freq
	if len(m) == 0 {
		return
	}
	var k int
	for pair := range m {
		k = len(pair.A)
		break
	}
	j := k - 1
	pd := newPairDeBruijn()
	for pair, mult := range m {
		if len(pair.A) != k || len(pair.B) != k {
			err = errors.New("kmers have different length")
			return
		}
		fr := pd.node(KmerPair{pair.A[:j], pair.B[:j]})
		to := pd.node(KmerPair{pair.A[1:], pair.B[1:]})
		for i := 0; i < mult; i++ {
			pd.arc(fr, to)
		}
	}
	return pd.g, ReadPairList{freq.D + 1, pd.jpairs}, nil
}

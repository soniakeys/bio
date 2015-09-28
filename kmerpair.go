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

func (pairs ReadPairList) OverlapSeq(order []int) (Seq, error) {
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
	jNode  map[KmerPair]int
	graph  [][]int
}

func newPairDeBruijn() *pairDeBruijn {
	return &pairDeBruijn{jNode: map[KmerPair]int{}}
}

func (pd *pairDeBruijn) node(jpair KmerPair) int {
	if n, ok := pd.jNode[jpair]; ok {
		return n
	}
	n := len(pd.jpairs)
	pd.jNode[jpair] = n
	pd.jpairs = append(pd.jpairs, jpair)
	pd.graph = append(pd.graph, nil)
	return n
}
func (pd *pairDeBruijn) arc(fr, to int) {
	pd.graph[fr] = append(pd.graph[fr], to)
}

// DeBruijn constructs a DeBruijn graph from a list of kmer pairs.
//
// The graph represents overlaps of length k-1 of the kmers.
// Nodes of the returned graph are labeled by
// k-1-mer (jmer) substring pairs of the kmer pairs.
func (pairs ReadPairList) DeBruijn() (graph [][]int, jpairs ReadPairList, err error) {
	ps := pairs.Pairs
	if len(ps) == 0 {
		return
	}
	k := len(ps[0].A)
	j := k - 1
	pd := newPairDeBruijn()
	for _, pair := range ps {
		if len(pair.A) != k || len(pair.B) != k {
			return nil, ReadPairList{},
				errors.New("kmers have different length")
		}
		fr := pd.node(KmerPair{pair.A[:j], pair.B[:j]})
		to := pd.node(KmerPair{pair.A[1:], pair.B[1:]})
		pd.arc(fr, to)
	}
	return pd.graph, ReadPairList{pairs.D + 1, pd.jpairs}, nil
}

func (pairs ReadPairList) Contigs() (cs []Seq, err error) {
	var g graph.AdjacencyList
	g, jpairs, err := pairs.DeBruijn()
	if err != nil {
		return nil, err
	}
	ps := g.MaximalNonBranchingPaths()
	for _, p := range ps {
		if s, err := jpairs.OverlapSeq(p); err == nil {
			cs = append(cs, s)
		} else {
			fmt.Println(err)
		}
	}
	return cs, nil
}

package bio

import (
	"bytes"
	"sort"
)

type taggedSuffix struct {
	suffix      []byte
	sourceIndex int
}

type ts []taggedSuffix

func (s ts) Len() int      { return len(s) }
func (s ts) Swap(i, j int) { s[i], s[j] = s[j], s[i] }
func (s ts) Less(i, j int) bool {
	return bytes.Compare(s[i].suffix, s[j].suffix) < 0
}

// SuffixArray digests multiple input sequences to provide a longest common
// sub sequence method.
type SuffixArray struct {
	a ts  // suffix "array"
	n int // number of sequences added
}

// AddSeq adds a sequence.  It can be called any number of times.
func (sa *SuffixArray) AddSeq(seq []byte) {
	for i := range seq {
		sa.a = append(sa.a, taggedSuffix{seq[i:], sa.n})
	}
	sa.n++
}

// LongestCommonSubSeq returns the longest common sub sequence of all added
// sequences.
func (sa *SuffixArray) LongestCommonSubSeq() []byte {
	switch sa.n {
	case 0:
		return nil
	case 1:
		return sa.a[0].suffix
	}
	// find k-longest common prefix of suffixes
	sort.Sort(sa.a)
	// populate sax with first suffix index of each source
	sax := make([]int, sa.n)
	for i := range sax {
		sax[i] = -1
	}
	for x, nFound := 0, 0; nFound < sa.n; x++ {
		xx := sa.a[x].sourceIndex
		if sax[xx] < 0 {
			sax[xx] = x
			nFound++
		}
	}
	// find common prefix of the first of each.
	// this becomes the initial longest common substring.
	lcs := lcp(sax, sa)
	// now search sa for longer ones, update s0 and preLen as they
	// are found
lcs:
	for x, s := range sa.a {
		// find next suffix to try replacing
		if x < sax[s.sourceIndex] {
			continue
		}
		// now search forward to find next suffix from this source
		y := x + 1
		for {
			if y == len(sa.a) {
				break lcs
			}
			if sa.a[y].sourceIndex == s.sourceIndex {
				break
			}
			y++
		}
		sax[s.sourceIndex] = y
		cp := lcp(sax, sa)
		if len(cp) > len(lcs) {
			lcs = cp
		}
	}
	return lcs
}

func lcp(sax []int, sa *SuffixArray) []byte {
	preLen := 0
	s0 := sa.a[sax[0]].suffix
lcp1:
	for preLen < len(s0) {
		for x := 1; x < len(sax); x++ {
			sx := sa.a[sax[x]].suffix
			if len(sx) == preLen {
				break lcp1
			}
			if sx[preLen] != s0[preLen] {
				break lcp1
			}
		}
		preLen++
	}
	return s0[:preLen]
}

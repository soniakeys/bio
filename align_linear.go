package bio

import (
	"bytes"
	"fmt"
)

// WIP
func AlignGlobalLS(s1, s2 Seq, a Aligner, indelPenalty int) (score int, a1, a2 Seq) {
	indel := -indelPenalty
	var lsa func(s1, s2 Seq)
	lsa = func(s1, s2 Seq) {
		if len(s1) == 0 {
			score += indel * len(s2)
			a1 = append(a1, bytes.Repeat([]byte{GapSymbol}, len(s2))...)
			a2 = append(a2, s2...)
			return
		}
		if len(s2) == 0 {
			score += indel * len(s1)
			a1 = append(a1, s1...)
			a2 = append(a2, bytes.Repeat([]byte{GapSymbol}, len(s1))...)
			return
		}
		p1, p2 := middleEdge(s1, s2, a, indelPenalty)
		lsa(s1[:p1.r], s2[:p1.c])
		switch {
		case p1.r == p2.r:
			score += indel
			a1 = append(a1, GapSymbol)
			a2 = append(a2, s2[p1.c])
		case p1.c == p2.c:
			score += indel
			a1 = append(a1, s1[p1.r])
			a2 = append(a2, GapSymbol)
		default:
			score += a.Score(s1[p1.r], s2[p1.c])
			a1 = append(a1, s1[p1.r])
			a2 = append(a2, s2[p1.c])
		}
		lsa(s1[p2.r:], s2[p2.c:])
	}
	lsa(s1, s2)
	return
}

type mPair struct{ r, c int }

func (p mPair) String() string { return fmt.Sprintf("(%d, %d)", p.r, p.c) }

func middleEdge(s1, s2 Seq, a Aligner, indelPenalty int) (p1, p2 mPair) {
	// some boundary cases
	if len(s1) == 0 {
		if len(s2) == 0 {
			panic("no middle edge with two null strings")
		}
		return mPair{0, 0}, mPair{0, 1}
	}
	if len(s2) == 0 {
		return mPair{0, 0}, mPair{1, 0}
	}
	indel := -indelPenalty
	mn, rev, rev1 := middleNodes(s1, s2, a, indel)
	//	fmt.Println("mn:", mn)
	last := len(rev) - 1
	max := 0
	set := false
	for _, p := range mn {
		//		fmt.Println("mn:", p)
		// gap s1
		g1 := rev1[last-p.r] + indel
		if !set || g1 > max {
			max = g1
			p1 = p
			p2 = mPair{p.r, p.c + 1}
		}
		g2 := rev[last-(p.r+1)] + indel
		if g2 > max {
			max = g2
			p1 = p
			p2 = mPair{p.r + 1, p.c}
		}
		//		ps := rev1[last-(p.r+1)]
		//		fmt.Println("ps:", ps)
		//		es := a.Score(s1[p.r], s2[p.c])
		//		fmt.Println("es:", es)
		mm := rev1[last-(p.r+1)] + a.Score(s1[p.r], s2[p.c])
		if mm > max {
			max = mm
			p1 = p
			p2 = mPair{p.r + 1, p.c + 1}
		}
	}
	return
}

func middleNodes(s1, s2 Seq, a Aligner, indel int) ([]mPair, []int, []int) {
	r1 := s1.Reverse()
	r2 := s2.Reverse()
	mc := len(s2) / 2
	fwd, _ := lastColumnScores(s1, s2[:mc], a, indel)
	rc := len(s2) - mc
	rev, rev1 := lastColumnScores(r1, r2[:rc], a, indel)
	//	fmt.Println("fwd:", fwd)
	//	fmt.Println("rev:", rev)
	last := len(rev) - 1
	r := []mPair{{0, mc}}
	max := fwd[0] + rev[last]
	for i := 1; i < len(fwd); i++ {
		switch s := fwd[i] + rev[last-i]; {
		case s > max:
			max = s
			r = []mPair{{i, mc}}
		case s == max:
			r = append(r, mPair{i, mc})
		}
	}
	return r, rev, rev1
}

// result sc0 is scores in last column, sc1 is scores in next to last clolumn
func lastColumnScores(s1, s2 Seq, a Aligner, indel int) (sc0, sc1 []int) {
	// score, base column.  a "column" is a slice parallel to s1.
	sc0 = make([]int, len(s1)+1)
	for s1x := range sc0 {
		sc0[s1x] = s1x * indel // s2 gaps, edge scores "down left side"
	}
	sc1 = make([]int, len(sc0))
	// propagate scores to last column
	for s2x := range s2 {
		sc1[0] = sc0[0] + indel // s1 gap "along top"
		for s1x := range s1 {
			max := sc1[s1x] + indel // s2 gap "down"
			if mm := sc0[s1x] + a.Score(s1[s1x], s2[s2x]); mm > max {
				max = mm // match/mismatch "diagonal"
			}
			if g := sc0[s1x+1] + indel; g > max {
				max = g // s1 gap "across"
			}
			sc1[s1x+1] = max
		}
		sc0, sc1 = sc1, sc0 // shift, recycle
	}
	return sc0, sc1
}

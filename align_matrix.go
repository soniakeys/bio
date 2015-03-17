package bio

import (
	"bytes"
	"fmt"
)

// Aligner defines a match-mismatch score for some alphabet.
type Aligner interface {
	Score(x, y byte) int
}

func AlignGlobal2(s1, s2 []byte, a Aligner, indelPenalty int) (score int, a1, a2 []byte) {
	indel := -indelPenalty
	var lsa func(s1, s2 []byte)
	lsa = func(s1, s2 []byte) {
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
		p1, p2 := MiddleEdge(s1, s2, a, indelPenalty)
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

type Pair struct{ r, c int }

func (p Pair) String() string { return fmt.Sprintf("(%d, %d)", p.r, p.c) }

func MiddleEdge(s1, s2 []byte, a Aligner, indelPenalty int) (p1, p2 Pair) {
	// some boundary cases
	if len(s1) == 0 {
		if len(s2) == 0 {
			panic("no middle edge with two null strings")
		}
		return Pair{0, 0}, Pair{0, 1}
	}
	if len(s2) == 0 {
		return Pair{0, 0}, Pair{1, 0}
	}
	indel := -indelPenalty
	mn, rev, rev1 := MiddleNodes(s1, s2, a, indel)
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
			p2 = Pair{p.r, p.c + 1}
		}
		g2 := rev[last-(p.r+1)] + indel
		if g2 > max {
			max = g2
			p1 = p
			p2 = Pair{p.r + 1, p.c}
		}
		//		ps := rev1[last-(p.r+1)]
		//		fmt.Println("ps:", ps)
		//		es := a.Score(s1[p.r], s2[p.c])
		//		fmt.Println("es:", es)
		mm := rev1[last-(p.r+1)] + a.Score(s1[p.r], s2[p.c])
		if mm > max {
			max = mm
			p1 = p
			p2 = Pair{p.r + 1, p.c + 1}
		}
	}
	return
}

func rev(s []byte) []byte {
	r := make([]byte, len(s))
	last := len(s) - 1
	for i, b := range s {
		r[last-i] = b
	}
	return r
}

func MiddleNodes(s1, s2 []byte, a Aligner, indel int) ([]Pair, []int, []int) {
	r1 := rev(s1)
	r2 := rev(s2)
	mc := len(s2) / 2
	fwd, _ := lastColumnScores(s1, s2[:mc], a, indel)
	rc := len(s2) - mc
	rev, rev1 := lastColumnScores(r1, r2[:rc], a, indel)
	//	fmt.Println("fwd:", fwd)
	//	fmt.Println("rev:", rev)
	last := len(rev) - 1
	r := []Pair{{0, mc}}
	max := fwd[0] + rev[last]
	for i := 1; i < len(fwd); i++ {
		switch s := fwd[i] + rev[last-i]; {
		case s > max:
			max = s
			r = []Pair{{i, mc}}
		case s == max:
			r = append(r, Pair{i, mc})
		}
	}
	return r, rev, rev1
}

// result sc0 is scores in last column, sc1 is scores in next to last clolumn
func lastColumnScores(s1, s2 []byte, a Aligner, indel int) (sc0, sc1 []int) {
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

// Align does global alignment using a dynamic programing technique.
func (m *SubstMatrix) Align(s1, s2 AA20, indelPenalty int) (a1, a2 AA, score int) {
	fmt.Println("Align", s1, s2, indelPenalty)
	s := make([][]int, len(s1)+1)
	b := make([][]byte, len(s1)+1)
	np := 0
	for i := range s {
		si := make([]int, len(s2)+1)
		si[0] = np
		np -= indelPenalty
		s[i] = si
		b[i] = make([]byte, len(s2)+1)
	}
	s0 := s[0]
	np = 0
	for j := range s0 {
		s0[j] = np
		np -= indelPenalty
	}
	var (
		// bb: backtrack code
		// '1' to match a symbol of s1
		// '2' to match a symbol of s2
		// '3' to match both
		bb  byte
		max int
	)
	for i := 1; i <= len(s1); i++ {
		for j := 1; j <= len(s2); j++ {
			fmt.Println("i, j:", i, j)
			for i, si := range s {
				fmt.Printf("%4d   %2d\n", si, b[i])
			}
			me := m.Score(s1[i-1], s2[j-1]) // me: matrix element
			sc3 := s[i-1][j-1] + me
			sc2 := s[i][j-1] - indelPenalty
			sc1 := s[i-1][j] - indelPenalty
			max = sc3
			bb = 3
			if sc2 > max {
				max = sc2
				bb = 2
			}
			if sc1 > max {
				max = sc1
				bb = 1
			}
			s[i][j] = max
			b[i][j] = bb
		}
	}
	score = s[len(s1)][len(s2)] // assign return value
	// backtrack to generate a1, a2 return values
	var c1, c2 byte
	var bt func(i, j, l int)
	bt = func(i, j, l int) {
		switch {
		case i == 0:
			a1 = make(AA, l+j)
			a2 = make(AA, l+j)
			for ; j > 0; j-- {
				a1[len(a1)-j] = GapSymbol
				a2[len(a2)-j] = s2[len(s2)-j]
			}
			return
		case j == 0:
			a1 = make(AA, l+i)
			a2 = make(AA, l+i)
			for ; i > 0; i-- {
				a1[len(a1)-i] = GapSymbol
				a2[len(a2)-i] = s1[len(s1)-i]
			}
			return
		}
		switch b[i][j] {
		case 1:
			bt(i-1, j, l+1)
			c1 = s1[i-1]
			c2 = GapSymbol
		case 2:
			bt(i, j-1, l+1)
			c1 = GapSymbol
			c2 = s2[j-1]
		default:
			bt(i-1, j-1, l+1)
			c1 = s1[i-1]
			c2 = s2[j-1]
		}
		a1[len(a1)-1-l] = c1
		a2[len(a2)-1-l] = c2
	}
	bt(len(s1), len(s2), 0)
	return
}

// LinearGap recomputes an alignment score.
//
// Arguments s and t are global alignments.  They must be of the same length.
//
// Result is score using linear gap pentalty gp.
func (m *SubstMatrix) LinearGap(s, t AA, gp int) (score int) {
	if len(s) != len(t) {
		panic("SubstMatrix.LinearGap strings different lengths")
	}
	for i, sa := range s {
		switch ta := t[i]; byte(GapSymbol) {
		case sa, ta:
			score -= gp
		default:
			score += m.Score(sa, ta)
		}
	}
	return
}

// ConstantGap recomputes an alignment score.
//
// Arguments s and t are global alignments.  They must be of the same length.
//
// Result is score using constant gap pentalty gp.
func (m *SubstMatrix) ConstantGap(s, t AA, gp int) (score int) {
	if len(s) != len(t) {
		panic("SubstMatrix.LinearGap strings different lengths")
	}
	const (
		noGap = iota
		sGap
		tGap
	)
	var openGap int
	for i, sa := range s {
		switch ta := t[i]; byte(GapSymbol) {
		case sa:
			if openGap != sGap {
				score -= gp
				openGap = sGap
			}
		case ta:
			if openGap != tGap {
				score -= gp
				openGap = tGap
			}
		default:
			score += m.Score(sa, ta)
			openGap = noGap
		}
	}
	return
}

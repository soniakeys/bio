package bio

var GapSymbol byte = '-' // represents a sequence alignment gap

// Aligner defines an integer match-mismatch score for some alphabet of bytes.
type Aligner interface {
	Score(x, y byte) int
}

// LinearGap recomputes an alignment score.
//
// Arguments s and t are alignment traces.  Gaps must be identified with
// the symbol of package variable `GapSymbol`.
//
// The function panics if s and t are not the same length.
//
// Result is score computed with linear gap pentalty.
func LinearGap(s, t Seq, a Aligner, gapPenalty int) (score int) {
	if len(s) != len(t) {
		panic("Sequences have different lengths")
	}
	for i, sa := range s {
		switch ta := t[i]; GapSymbol {
		case sa, ta:
			score -= gapPenalty
		default:
			score += a.Score(sa, ta)
		}
	}
	return
}

// ConstantGap recomputes an alignment score.
//
// Arguments s and t are alignment traces.  Gaps must be identified with
// the symbol of package variable `GapSymbol`.
//
// The function panics if s and t are not the same length.
//
// Result is score computed with constant gap pentalty.
func ConstantGap(s, t Seq, a Aligner, gapPenalty int) (score int) {
	if len(s) != len(t) {
		panic("Sequences have different lengths")
	}
	const (
		noGap = iota
		sGap
		tGap
	)
	var openGap int
	for i, sa := range s {
		switch ta := t[i]; GapSymbol {
		case sa:
			if openGap != sGap {
				score -= gapPenalty
				openGap = sGap
			}
		case ta:
			if openGap != tGap {
				score -= gapPenalty
				openGap = tGap
			}
		default:
			score += a.Score(sa, ta)
			openGap = noGap
		}
	}
	return
}

// AlignGlobal.
//
// Algorithm is Needleman–Wunsch.
//
// Result traces t1, t2 use the package variable GapSymbol to indicate gaps.
func AlignGlobal(s1, s2 Seq, a Aligner, indelPenalty int) (score int, t1, t2 Seq) {
	stride := len(s2) + 1
	// score matrix/node labels
	s := make([]int, (len(s1)+1)*stride)
	// backtrack matrix.  parallel to s.
	b := make([]int, len(s))
	// values to store in b representing gap in s1, gap in s2, and
	// match/mismatch respectively.  these can be subtracted from a
	// position/node number to get the previous position/node number.
	g1 := 1
	g2 := stride
	mm := stride + 1
	// index into s, b
	x := 0

	// "top row" gap in s1 all across
	for x = 1; x < stride; x++ {
		s[x] = s[x-g1] - indelPenalty
		b[x] = g1
	}
	// "lower rows"
	for i := range s1 {
		// first position is gap in s2
		s[x] = s[x-g2] - indelPenalty
		b[x] = g2
		x++
		// "interior" positions
		for j := range s2 {
			// case corresponds to gap in s2
			sMax := s[x-g2] - indelPenalty
			bMax := g2
			// case corresponds to gap in s1
			if s0 := s[x-g1] - indelPenalty; s0 > sMax {
				sMax = s0
				bMax = g1
			}
			// case corresponds to match/mismatch
			if s0 := s[x-mm] + a.Score(s1[i], s2[j]); s0 > sMax {
				sMax = s0
				bMax = mm
			}
			// store accumulated max of three cases
			s[x] = sMax
			b[x] = bMax
			x++
		}
	}
	x = len(s) - 1
	score = s[x]
	// backtrack then reverse
	for ; x > 0; x -= b[x] {
		switch b[x] {
		case g1:
			t1 = append(t1, GapSymbol)
			t2 = append(t2, s2[x%stride-1])
		case g2:
			t1 = append(t1, s1[x/stride-1])
			t2 = append(t2, GapSymbol)
		default:
			t1 = append(t1, s1[x/stride-1])
			t2 = append(t2, s2[x%stride-1])
		}
	}
	last := len(t1) - 1
	for i := range t1[:len(t1)/2] {
		t1[i], t1[last-i] = t1[last-i], t1[i]
		t2[i], t2[last-i] = t2[last-i], t2[i]
	}
	return
}

// AlignLocal.
//
// Algorithm is Smith–Waterman.
//
// Result traces t1, t2 use the package variable GapSymbol to indicate gaps.
func AlignLocal(s1, s2 Seq, a Aligner, indelPenalty int) (score int, t1, t2 Seq) {
	stride := len(s2) + 1
	// score matrix/node labels
	s := make([]int, (len(s1)+1)*stride)
	// backtrack matrix.  parallel to s.
	b := make([]byte, len(s))
	// values to store in b representing a choice made for a node
	const (
		sp byte = 'p' // skip prefix
		g1 byte = '1' // gap in s1
		g2 byte = '2' // gap in s2
		mm byte = 'm' // match/mismatch
	)
	// values to subtract from a position/node number to get the previous
	// position/node number.
	g1x := 1
	g2x := stride
	mmx := stride + 1

	// index into s, b
	x := 0

	// "top row" skip prefix all across
	for x = 1; x < stride; x++ {
		b[x] = sp
	}
	// "lower rows"
	for i := range s1 {
		// first position skips prefix
		b[x] = sp
		x++
		// "interior" positions
		for j := range s2 {
			// case corresponds to skipping the prefix, a "free taxi ride"
			sMax := 0
			bMax := sp
			// case corresponds to gap in s2
			if s0 := s[x-g2x] - indelPenalty; s0 > sMax {
				sMax = s0
				bMax = g2
			}
			// case corresponds to gap in s1
			if s0 := s[x-g1x] - indelPenalty; s0 > sMax {
				sMax = s0
				bMax = g1
			}
			// case corresponds to match/mismatch
			if s0 := s[x-mmx] + a.Score(s1[i], s2[j]); s0 > sMax {
				sMax = s0
				bMax = mm
			}
			// store accumulated max of three cases
			s[x] = sMax
			b[x] = bMax
			x++
		}
	}
	sLast := len(s) - 1
	score = s[sLast] // final score
	// well, except there's a different rule for the last node.
	// preference goes to skipping a suffix.  so the max score out of
	// all of s can be picked as the last node.
	score = 0 // accumulate best score for result
	x = 0     // also accumulate corresponding node number
	for xx, sx := range s {
		if sx > score {
			score = sx
			x = xx
		}
	}
	// x becomes the starting point for backtracking
	// backtrack then reverse
bt:
	for x > 0 {
		switch b[x] {
		case g1:
			t1 = append(t1, GapSymbol)
			t2 = append(t2, s2[x%stride-1])
			x -= g1x
		case g2:
			t1 = append(t1, s1[x/stride-1])
			t2 = append(t2, GapSymbol)
			x -= g2x
		case mm:
			t1 = append(t1, s1[x/stride-1])
			t2 = append(t2, s2[x%stride-1])
			x -= mmx
		default:
			break bt
		}
	}
	last := len(t1) - 1
	for i := range t1[:len(t1)/2] {
		t1[i], t1[last-i] = t1[last-i], t1[i]
		t2[i], t2[last-i] = t2[last-i], t2[i]
	}
	return
}

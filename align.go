package bio

func AlignGlobal(s1, s2 Seq, a Aligner, indelPenalty int) (score int, t1, t2 Seq) {
	s11 := len(s1) + 1
	s21 := len(s2) + 1
	// score matrix/node labels
	s := make([]int, s11*s21)
	// backtrack matrix.  parallel to s.
	b := make([]int, len(s))
	// values to store in b representing gap in s1, gap in s2, and
	// match/mismatch respectively.  these can be subtracted from a
	// position/node number to get the previous position/node number.
	g1 := 1
	g2 := s21
	mm := s21 + 1
	// index into s, b
	x := 0

	// "top row" gap in s1 all across
	for x = 1; x < s21; x++ {
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
			t1 = append(t1, '-')
			t2 = append(t2, s2[x%s21-1])
		case g2:
			t1 = append(t1, s1[x/s21-1])
			t2 = append(t2, '-')
		default:
			t1 = append(t1, s1[x/s21-1])
			t2 = append(t2, s2[x%s21-1])
		}
	}
	last := len(t1) - 1
	for i := range t1[:len(t1)/2] {
		t1[i], t1[last-i] = t1[last-i], t1[i]
		t2[i], t2[last-i] = t2[last-i], t2[i]
	}
	return
}

func AlignLocal(s1, s2 Seq, a Aligner, indelPenalty int) (score int, t1, t2 Seq) {
	s11 := len(s1) + 1
	s21 := len(s2) + 1
	// score matrix/node labels
	s := make([]int, s11*s21)
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
	g2x := s21
	mmx := s21 + 1

	// index into s, b
	x := 0

	// "top row" skip prefix all across
	for x = 1; x < s21; x++ {
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
	// preference goes to skipping a suffix.  to fix up the computation
	// done so far for the last node, here compute the best suffix
	// then use it if it's no worse than the score computed above.
	sMax := 0 // best score
	xMax := 0 // node number giving best score
	for x, sx := range s[:sLast] {
		if sx > sMax {
			sMax = sx
			xMax = x
		}
	}
	x = sLast // x becomes the starting point for backtracking
	if sMax >= score {
		score = sMax
		x = xMax
	}
	// backtrack then reverse
bt:
	for x > 0 {
		switch b[x] {
		case g1:
			t1 = append(t1, '-')
			t2 = append(t2, s2[x%s21-1])
			x -= g1x
		case g2:
			t1 = append(t1, s1[x/s21-1])
			t2 = append(t2, '-')
			x -= g2x
		case mm:
			t1 = append(t1, s1[x/s21-1])
			t2 = append(t2, s2[x%s21-1])
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

func AlignPair(mode string, s1, s2 Seq, a Aligner, indelPenalty int) (score int, t1, t2 Seq) {
	pa := newPairAligner(s1, s2, a, indelPenalty)
	switch mode {
	case "global":
		pa.setGlobal()
	case "local":
		pa.setLocal()
	default:
		return -1, nil, nil
	}
	pa.align()
	score = pa.score
	t1, t2 = pa.trace()
	return
}

type btFunc func(x int) (px int, b1, b2 byte)

type pairAligner struct {
	// parameters
	s1, s2       Seq
	a            Aligner
	indelPenalty int

	// a few values trivially derived from parameters but assigned to
	// identifier names more meaningful in certain contexts.
	s21 int // = len(s2) + 1, stride for s, b.
	g1x int // x increment for s1 gap, = 1
	g2x int // x increment for s2 gap, = s21
	mmx int // x increment for match/mismatch, = s21 + 1

	// configuration for alignment mode
	top, left *paRule
	interior  []*paRule
	final     func()

	// main dynamic programming representation
	s []int    // score matrix/nodes
	b []btFunc // backtrack matrix.  parallel to s.

	// results
	xLast int // alignment end, backtrack start. used by trace()
	// score is final alignment score.
	// should be s[len(s)-1] but this is clearer.
	score int
}

func newPairAligner(s1, s2 Seq, a Aligner, indelPenalty int) *pairAligner {
	pa := &pairAligner{s1: s1, s2: s2, a: a, indelPenalty: indelPenalty}
	pa.s = make([]int, (len(s1)+1)*(len(s2)+1))
	pa.b = make([]btFunc, len(pa.s))
	pa.s21 = len(s2) + 1
	pa.g1x = 1
	pa.g2x = pa.s21
	pa.mmx = pa.s21 + 1
	return pa
}

// skip prefix rule, "free taxi ride" (free taxi from start, that is)
func (pa *pairAligner) spScore(x int) int              { return 0 }
func (pa *pairAligner) spBack(x int) (int, byte, byte) { return 0, 0, 0 }

// gap in s2
func (pa *pairAligner) g2Score(x int) int {
	return pa.s[x-pa.g2x] - pa.indelPenalty
}
func (pa *pairAligner) g2Back(x int) (int, byte, byte) {
	return x - pa.g2x, pa.s1[x/pa.s21-1], '-'
}

// gap in s1
func (pa *pairAligner) g1Score(x int) int {
	return pa.s[x-pa.g1x] - pa.indelPenalty
}
func (pa *pairAligner) g1Back(x int) (int, byte, byte) {
	return x - pa.g1x, '-', pa.s2[x%pa.s21-1]
}

// match/mismatch
func (pa *pairAligner) mmScore(x int) int {
	i := (x / pa.s21) - 1
	j := (x % pa.s21) - 1
	return pa.s[x-pa.mmx] + pa.a.Score(pa.s1[i], pa.s2[j])
}
func (pa *pairAligner) mmBack(x int) (int, byte, byte) {
	return x - pa.mmx, pa.s1[x/pa.s21-1], pa.s2[x%pa.s21-1]
}

// initialize with method values of above method
type paRule struct {
	score func(x int) int // score for indicated choice
	// back one step from x.
	// returns previous x, symbols to add to traces t1, t2
	back func(x int) (px int, b1, b2 byte)
}

func (pa *pairAligner) setGlobal() {
	// "top row", all gaps in s1
	pa.top = &paRule{pa.g1Score, pa.g1Back}
	// "left edge", all gaps in s2
	pa.left = &paRule{pa.g2Score, pa.g2Back}
	// order of rules important for repeatability.
	// a different order may give different traces, although with
	// the same score.
	pa.interior = []*paRule{
		&paRule{pa.g2Score, pa.g2Back}, // gap in s2
		&paRule{pa.g1Score, pa.g1Back}, // gap in s1
		&paRule{pa.mmScore, pa.mmBack}, // match/mismatch
	}
	// no extra work to do except just set final result values
	pa.final = func() {
		pa.xLast = len(pa.s) - 1
		pa.score = pa.s[pa.xLast]
	}
}

func (pa *pairAligner) setLocal() {
	pa.top = &paRule{pa.spScore, pa.spBack}
	pa.left = &paRule{pa.spScore, pa.spBack}
	pa.interior = []*paRule{
		&paRule{pa.spScore, pa.spBack}, // skip prefix
		&paRule{pa.g2Score, pa.g2Back}, // gap in s2
		&paRule{pa.g1Score, pa.g1Back}, // gap in s1
		&paRule{pa.mmScore, pa.mmBack}, // match/mismatch
	}
	pa.final = func() {
		sLast := len(pa.s) - 1
		// there's a different rule for the last node.
		// preference goes to skipping a suffix.  to fix up the computation
		// done so far for the last node, here compute the best suffix
		// then use it if it's no worse than the score computed above.
		sMax := 0 // best score
		xMax := 0 // node number giving best score
		for x, sx := range pa.s[:sLast] {
			if sx > sMax {
				sMax = sx
				xMax = x
			}
		}
		pa.xLast = sLast
		if sMax >= pa.s[sLast] {
			pa.s[sLast] = sMax
			pa.xLast = xMax
		}
		pa.score = pa.s[sLast]
	}
}

func (pa *pairAligner) trace() (t1, t2 Seq) {
	// backtrack then reverse
	x := pa.xLast
	for x > 0 {
		px, b1, b2 := pa.b[x](x)
		if b1 != 0 {
			t1 = append(t1, b1)
		}
		if b2 != 0 {
			t2 = append(t2, b2)
		}
		x = px
	}
	last := len(t1) - 1
	for i := range t1[:len(t1)/2] {
		t1[i], t1[last-i] = t1[last-i], t1[i]
		t2[i], t2[last-i] = t2[last-i], t2[i]
	}
	return
}

func (pa *pairAligner) align() {
	x := 1 // node number.  index into s, b
	for range pa.s2 {
		pa.s[x] = pa.top.score(x)
		pa.b[x] = pa.top.back
		x++
	}
	// "lower rows"
	for range pa.s1 {
		pa.s[x] = pa.left.score(x)
		pa.b[x] = pa.left.back
		x++
		// "interior" positions
		for range pa.s2 {
			sMax := pa.interior[0].score(x)
			bMax := pa.interior[0].back
			for _, r := range pa.interior[1:] {
				if s0 := r.score(x); s0 > sMax {
					sMax = s0
					bMax = r.back
				}
			}
			pa.s[x] = sMax
			pa.b[x] = bMax
			x++
		}
	}
	pa.final()
}

package bio

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
			t1 = append(t1, '-')
			t2 = append(t2, s2[x%stride-1])
		case g2:
			t1 = append(t1, s1[x/stride-1])
			t2 = append(t2, '-')
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
	for xx, sx := range s[1:] {
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
			t1 = append(t1, '-')
			t2 = append(t2, s2[x%stride-1])
			x -= g1x
		case g2:
			t1 = append(t1, s1[x/stride-1])
			t2 = append(t2, '-')
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

func AlignPair(mode string, s1, s2 Seq, a Aligner, indelPenalty int) (score int, t1, t2 Seq) {
	pa := newPairAligner(s1, s2, a, indelPenalty)
	switch mode {
	case "global":
		pa.setGlobal()
	case "local":
		pa.setLocal()
	case "fitting":
		pa.setFitting()
	case "overlap":
		pa.setOverlap()
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

	// a couple values trivially derived from parameters but assigned to
	// identifier names more meaningful in certain contexts.
	stride int // = len(s2) + 1, stride for s, b.
	mmx    int // x increment for match/mismatch, = stride + 1

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
	pa.stride = len(s2) + 1
	pa.mmx = pa.stride + 1 // match/mismatch increment
	return pa
}

// skip prefix rule, "free taxi ride" (free taxi from start, that is)
func spScore(x int) int                 { return 0 }
func spBack(x int) (int, byte, byte)    { return 0, 0, 0 }
func (pa *pairAligner) spRule() *paRule { return &paRule{spScore, spBack} }

// gap in s2
func (pa *pairAligner) g2Score(x int) int {
	return pa.s[x-pa.stride] - pa.indelPenalty
}
func (pa *pairAligner) g2Back(x int) (int, byte, byte) {
	return x - pa.stride, pa.s1[x/pa.stride-1], '-'
}
func (pa *pairAligner) g2Rule() *paRule {
	return &paRule{pa.g2Score, pa.g2Back}
}

// gap in s1
func (pa *pairAligner) g1Score(x int) int {
	return pa.s[x-1] - pa.indelPenalty
}
func (pa *pairAligner) g1Back(x int) (int, byte, byte) {
	return x - 1, '-', pa.s2[x%pa.stride-1]
}
func (pa *pairAligner) g1Rule() *paRule {
	return &paRule{pa.g1Score, pa.g1Back}
}

// match/mismatch
func (pa *pairAligner) mmScore(x int) int {
	i := (x / pa.stride) - 1
	j := (x % pa.stride) - 1
	return pa.s[x-pa.mmx] + pa.a.Score(pa.s1[i], pa.s2[j])
}
func (pa *pairAligner) mmBack(x int) (int, byte, byte) {
	return x - pa.mmx, pa.s1[x/pa.stride-1], pa.s2[x%pa.stride-1]
}
func (pa *pairAligner) mmRule() *paRule {
	return &paRule{pa.mmScore, pa.mmBack}
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
	pa.top = pa.g1Rule()
	// "left edge", all gaps in s2
	pa.left = pa.g2Rule()
	// order of rules can matter.
	// a different order may give different traces, although with
	// the same score.
	pa.interior = []*paRule{
		pa.g2Rule(), // gap in s2
		pa.g1Rule(), // gap in s1
		pa.mmRule(), // match/mismatch
	}
	// no extra work to do except just set final result values
	pa.final = func() {
		pa.xLast = len(pa.s) - 1
		pa.score = pa.s[pa.xLast]
	}
}

func (pa *pairAligner) setLocal() {
	pa.top = pa.spRule()  // skip prefix
	pa.left = pa.spRule() // skip prefix
	pa.interior = []*paRule{
		pa.spRule(), // skip prefix
		pa.g2Rule(), // gap in s2
		pa.g1Rule(), // gap in s1
		pa.mmRule(), // match/mismatch
	}
	pa.final = func() {
		// rule for the last node:
		// preference goes to skipping a suffix.
		// max score out of all of s becomes the last node.
		sMax := 0 // best score
		xMax := 0 // node number giving best score
		for x, sx := range pa.s[1:] {
			if sx > sMax {
				sMax = sx
				xMax = x
			}
		}
		/* alternatively, by "columns"
		for j := 0; j < pa.stride; j++ {
			for x := j; x < len(pa.s); x += pa.stride {
				if sx := pa.s[x]; sx > sMax {
					sMax = sx
					xMax = x
				}
			}
		}
		*/
		pa.xLast = xMax
		pa.score = sMax
	}
}

// fit s1 (shorter seq) within s2 (longer seq)
// so "matrix" is wider than it is tall
func (pa *pairAligner) setFitting() {
	pa.top = pa.spRule()  // sp, top row like local
	pa.left = pa.g2Rule() // g2, left edge like global
	pa.interior = []*paRule{
		pa.g2Rule(), // gap in s2
		pa.g1Rule(), // gap in s1
		pa.mmRule(), // match/mismatch
	}
	pa.final = pa.skipS2Suffix
}

// a little method shared by setFitting and setOverlap
func (pa *pairAligner) skipS2Suffix() {
	// skip a suffix of s2 only.
	// fix up is to take the best score "across the bottom"
	xMax := len(pa.s) - pa.stride // node number giving best score
	sMax := pa.s[xMax]            // best score
	for x := xMax + 1; x < len(pa.s); x++ {
		if sx := pa.s[x]; sx > sMax {
			sMax = sx
			xMax = x
		}
	}
	pa.xLast = xMax
	pa.score = sMax
}

// align a suffix of s1 with a prefix of s2.
// in other words, skip a prefix of s1, skip a suffix of s2.
func (pa *pairAligner) setOverlap() {
	pa.top = pa.g1Rule()  // g1, top row like global
	pa.left = pa.spRule() // sp, left edge like local
	pa.interior = []*paRule{
		pa.g2Rule(), // gap in s2
		pa.g1Rule(), // gap in s1
		pa.mmRule(), // match/mismatch
	}
	pa.final = pa.skipS2Suffix
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

package bio

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
	return x - pa.stride, pa.s1[x/pa.stride-1], GapSymbol
}
func (pa *pairAligner) g2Rule() *paRule {
	return &paRule{pa.g2Score, pa.g2Back}
}

// gap in s1
func (pa *pairAligner) g1Score(x int) int {
	return pa.s[x-1] - pa.indelPenalty
}
func (pa *pairAligner) g1Back(x int) (int, byte, byte) {
	return x - 1, GapSymbol, pa.s2[x%pa.stride-1]
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
		for x, sx := range pa.s {
			if sx > sMax {
				sMax = sx
				xMax = x
			}
		}
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

/*
func AlignAffine(mode string, s1, s2 Seq, a Aligner, gapOpenPenalty, gapExendPenalty float64) (score float64, t1, t2 Seq) {
	aa := newAffineAligner(s1, s2, a, gapOpenPenalty, gapExendPenalty)
	fmt.Println("AffineAlign len(s1), len(s2) =", len(s1), len(s2))
	switch mode {
	case "global":
		aa.setGlobal()
			case "local":
				aa.setLocal()
			case "fitting":
				aa.setFitting()
			case "overlap":
				aa.setOverlap()
	default:
		return -1, nil, nil
	}
	aa.align()
	score = aa.score
	t1, t2 = aa.trace()
	return
}

type affineAligner struct {
	// parameters
	s1, s2           Seq
	a                Aligner
	gapOpenPenalty   float64
	gapExtendPenalty float64

	// a couple values trivially derived from parameters but assigned to
	// identifier names more meaningful in certain contexts.
	stride int // = len(s2) + 1, stride for s, b.
	mmx    int // x increment for match/mismatch, = stride + 1

	// configuration for alignment mode
	g1Rules []*aaRule
	g2Rules []*aaRule
	mmRules []*aaRule
	final   func()

	// main dynamic programming representation
	// score matrix/nodes; for gap s1, gap s2, match/mismatch levels
	sg1, sg2, smm []float64
	bg1, bg2, bmm []abFunc // backtrack matrix.  parallel to s.

	// results
	xLast int     // node where alignment ends. used by trace
	fLast abFunc  // bt func where alignment ends. used by trace
	score float64 // score is final alignment score.
}

// back one step from x.
// returns previous x, abFunc from previous level at previous x,
// symbols to add to traces t1, t2
type abFunc func(x int) (px int, pf abFunc, b1, b2 byte, pLev string)

func newAffineAligner(s1, s2 Seq, a Aligner, gapOpenPenalty, gapExtendPenalty float64) *affineAligner {
	aa := &affineAligner{
		s1: s1, s2: s2, a: a,
		gapOpenPenalty:   gapOpenPenalty,
		gapExtendPenalty: gapExtendPenalty,
	}
	aa.sg1 = make([]float64, (len(s1)+1)*(len(s2)+1))
	aa.sg2 = make([]float64, len(aa.sg1))
	aa.smm = make([]float64, len(aa.sg1))
	aa.bg1 = make([]abFunc, len(aa.sg1))
	aa.bg2 = make([]abFunc, len(aa.sg1))
	aa.bmm = make([]abFunc, len(aa.sg1))
	aa.stride = len(s2) + 1
	aa.mmx = aa.stride + 1 // match/mismatch increment
	return aa
}

type aaRule struct {
	score func(x int) float64
	back  abFunc
	name  string // for debugging :(
}

// open gap in s2
func (aa *affineAligner) g2OpenScore(x int) float64 {
	px := x - aa.stride
	if px < 0 {
		return math.Inf(-1)
	}
	return aa.smm[px] - aa.gapOpenPenalty
}
func (aa *affineAligner) g2OpenBack(x int) (int, abFunc, byte, byte, string) {
	px := x - aa.stride
	return px, aa.bmm[px], aa.s1[px/aa.stride], GapSymbol, "mm"
}
func (aa *affineAligner) g2OpenRule() *aaRule {
	return &aaRule{aa.g2OpenScore, aa.g2OpenBack, "g2 open"}
}

// extend gap in s2
func (aa *affineAligner) g2ExtendScore(x int) float64 {
	px := x - aa.stride
	if px < aa.stride {
		return math.Inf(-1)
	}
	return aa.sg2[px] - aa.gapExtendPenalty
}
func (aa *affineAligner) g2ExtendBack(x int) (int, abFunc, byte, byte, string) {
	px := x - aa.stride
	return px, aa.bg2[px], aa.s1[px/aa.stride], GapSymbol, "g2"
}
func (aa *affineAligner) g2ExtendRule() *aaRule {
	return &aaRule{aa.g2ExtendScore, aa.g2ExtendBack, "g2 extend"}
}

// close gap in s2
func (aa *affineAligner) g2CloseScore(x int) float64 {
	return aa.sg2[x]
}
func (aa *affineAligner) g2CloseBack(x int) (int, abFunc, byte, byte, string) {
	return x, aa.bg2[x], 0, 0, "g2"
}
func (aa *affineAligner) g2CloseRule() *aaRule {
	return &aaRule{aa.g2CloseScore, aa.g2CloseBack, "g2 close"}
}

// open gap in s1
func (aa *affineAligner) g1OpenScore(x int) float64 {
	if cx := x % aa.stride; cx == 0 {
		return math.Inf(-1)
	}
	return aa.smm[x-1] - aa.gapOpenPenalty
}
func (aa *affineAligner) g1OpenBack(x int) (int, abFunc, byte, byte, string) {
	px := x - 1
	return px, aa.bmm[px], GapSymbol, aa.s2[px%aa.stride], "mm"
}
func (aa *affineAligner) g1OpenRule() *aaRule {
	return &aaRule{aa.g1OpenScore, aa.g1OpenBack, "g1 open"}
}

// extend gap in s1
func (aa *affineAligner) g1ExtendScore(x int) float64 {
	if cx := x % aa.stride; cx <= 1 {
		return math.Inf(-1)
	}
	return aa.sg1[x-1] - aa.gapExtendPenalty
}
func (aa *affineAligner) g1ExtendBack(x int) (int, abFunc, byte, byte, string) {
	px := x - 1
	return px, aa.bg1[px], GapSymbol, aa.s2[px%aa.stride], "g1"
}
func (aa *affineAligner) g1ExtendRule() *aaRule {
	return &aaRule{aa.g1ExtendScore, aa.g1ExtendBack, "g1 extend"}
}

// close gap in s1
func (aa *affineAligner) g1CloseScore(x int) float64 {
	return aa.sg1[x]
}
func (aa *affineAligner) g1CloseBack(x int) (int, abFunc, byte, byte, string) {
	return x, aa.bg1[x], 0, 0, "g1"
}
func (aa *affineAligner) g1CloseRule() *aaRule {
	return &aaRule{aa.g1CloseScore, aa.g1CloseBack, "g1 close"}
}

// match/mismatch
func (aa *affineAligner) mmScore(x int) float64 {
	i := (x / aa.stride) - 1
	j := (x % aa.stride) - 1
	if i < 0 || j < 0 {
		return math.Inf(-1)
	}
	fmt.Printf("      score %c %c = %d\n",
		aa.s1[i], aa.s2[j], aa.a.Score(aa.s1[i], aa.s2[j]))
	return aa.smm[x-aa.mmx] + float64(aa.a.Score(aa.s1[i], aa.s2[j]))
}
func (aa *affineAligner) mmBack(x int) (int, abFunc, byte, byte, string) {
	px := x - aa.mmx
	return px, aa.bmm[x], aa.s1[px/aa.stride], aa.s2[px%aa.stride], "mm"
}
func (aa *affineAligner) mmRule() *aaRule {
	return &aaRule{aa.mmScore, aa.mmBack, "match/mismatch"}
}

func (aa *affineAligner) setGlobal() {
	fmt.Println("setGlobal aa.stride", aa.stride)
	aa.g1Rules = []*aaRule{
		aa.g1ExtendRule(), // extend previous gap
		aa.g1OpenRule(),   // open new one
	}
	aa.g2Rules = []*aaRule{
		aa.g2ExtendRule(),
		aa.g2OpenRule(),
	}
	aa.mmRules = []*aaRule{
		aa.g1CloseRule(),
		aa.mmRule(),
		aa.g2CloseRule(),
	}
	// final result is level with max score in last position
	aa.final = func() {
		last := len(aa.sg1) - 1
		sMax := aa.sg1[last]
		fMax := aa.bg1[last]
		if s := aa.smm[last]; s > sMax {
			sMax = s
			fMax = aa.bmm[last]
		}
		if s := aa.sg2[last]; s > sMax {
			sMax = s
			fMax = aa.bg2[last]
		}
		aa.score = sMax
		aa.fLast = fMax
		aa.xLast = last
	}
}

func applyRules(rules []*aaRule, x int) (float64, abFunc) {
	fmt.Println("apply rules node", x, len(rules), "rules")
	sMax := rules[0].score(x)
	bMax := rules[0].back
	nMax := rules[0].name
	fmt.Println("  ", nMax, "score", sMax)
	for _, r := range rules[1:] {
		s0 := r.score(x)
		fmt.Println("  ", r.name, "score", s0)
		if s0 > sMax {
			sMax = s0
			bMax = r.back
			nMax = r.name
		}
	}
	fmt.Println("max is", nMax, "score", sMax)
	return sMax, bMax
}

func (aa *affineAligner) align() {
	for x := 1; x < len(aa.sg1); x++ {
		aa.sg1[x], aa.bg1[x] = applyRules(aa.g1Rules, x)
		aa.sg2[x], aa.bg2[x] = applyRules(aa.g2Rules, x)
		aa.smm[x], aa.bmm[x] = applyRules(aa.mmRules, x)
	}
	aa.final()
	aa.dump()
}

func (aa *affineAligner) trace() (t1, t2 Seq) {
	// backtrack then reverse
	x := aa.xLast
	f := aa.fLast
	fmt.Println("trace! xLast, fLast", x, f)
	var b1, b2 byte
	for x > 0 {
		fmt.Printf("x, f: %d %x ", x, f)
		x, f, b1, b2, _ = f(x)
		fmt.Printf("-> b1, b2: %c %c\n", b1, b2)
		if b1 != 0 {
			t1 = append(t1, b1)
			t2 = append(t2, b2)
		}
	}
	last := len(t1) - 1
	for i := range t1[:len(t1)/2] {
		t1[i], t1[last-i] = t1[last-i], t1[i]
		t2[i], t2[last-i] = t2[last-i], t2[i]
	}
	return
}

func (aa *affineAligner) dump() {
	fmt.Println("dump:")
	fmt.Print("g1   ")
	for _, b := range aa.s2 {
		fmt.Print("  ", string(b))
	}
	for x, _ := range aa.sg1 {
		// s1 symbol in left column
		if x%aa.stride == 0 {
			fmt.Println()
			if x == 0 {
				fmt.Print(" ")
			} else {
				i := x / aa.stride
				fmt.Print(string(aa.s1[i-1]))
			}
		}
		// row 1 of 3: score, bt dir.
		pLev := "--"
		if bf := aa.bg1[x]; bf != nil {
			_, _, _, _, pLev = bf(x)
		}
		fmt.Printf("%6.0f %s", aa.sg1[x], pLev)
	}
}
*/

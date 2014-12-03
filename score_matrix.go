package bio

// ScoreMatrix holds an score matrix used for protein sequence alignment.
type ScoreMatrix struct {
	posTable   [25]int
	scoreTable [20 * 20]int
}

// NewScoreMatrix creates a new ScoreMatrix given alphabet a and score table s
// as a flat slice.
//
// The alphabet a must contain exactly the symbols of AA20Alphabet but
// in the order of the axes of the score table.
// See blosum62.go and pam250.go for examples.
//
// Unless otherwise documented, the allowable alphabet for AA strings passed
// to methods is AA20Alphabet + GapSymbol.  Other symbols may cause panic.
func NewScoreMatrix(a string, s []int) *ScoreMatrix {
	if len(a) != 20 || len(s) != 400 {
		return nil
	}
	m := &ScoreMatrix{}
	for i, aa := range a {
		m.posTable[aa-'A'] = i
	}
	copy(m.scoreTable[:], s)
	return m
}

// returns row position in score table
func (m *ScoreMatrix) pos(aa byte) int {
	return m.posTable[aa-'A']
}

// Score returns the alignment score of the given amino acid symbols.
//
// Arguments must be in AA20Alphabet.  GapSymbol will cause panic.
func (m *ScoreMatrix) Score(a1, a2 byte) int {
	return m.scoreTable[m.pos(a1)*20+m.pos(a2)]
}

func (m *ScoreMatrix) LinearGap(s, t AA, gp int) (score int) {
	if len(s) != len(t) {
		panic("ScoreMatrix.LinearGap strings different lengths")
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

func (m *ScoreMatrix) ConstantGap(s, t AA, gp int) (score int) {
	if len(s) != len(t) {
		panic("ScoreMatrix.LinearGap strings different lengths")
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

package bio

import (
	"fmt"
	"math"
)

// s1 is v in the text, indexed by i, and labeling down the left of the matrix
// s2 is w in the text, indexed by j, and labeling across the top

// WIP
func AlignGlobalAffine(s1, s2 Seq, a Aligner, gapOpenPenalty, gapExtendPenalty float64) (score float64, t1, t2 Seq) {
	stride := len(s2) + 1
	// score matrices
	// sg1 holds scores for the state where gap in s1 is open
	sg1 := make([]float64, (len(s1)+1)*stride) // "upper level" by the text
	smm := make([]float64, (len(s1)+1)*stride) // "middle level"
	sg2 := make([]float64, (len(s1)+1)*stride) // "lower level"
	// values to store in b representing gap in s1, gap in s2, and
	// match/mismatch respectively.  these can be subtracted from a
	// position/node number to get the previous position/node number.
	const g1 = 1
	g2 := stride // g2, mm are also constant after set
	mm := stride + 1
	// backtrack info is a level-index pair
	type bi struct {
		level int // previous level. use g1, g2, mm values just defined
		x     int // previous x.
	}
	// backtrack matrices.  parallel to s.
	bg1 := make([]bi, len(sg1))
	bg2 := make([]bi, len(sg1))
	bmm := make([]bi, len(sg1))

	// start of alignment has no gaps open, so -inf score for gapped states
	sg1[0] = math.Inf(-1)
	sg2[0] = math.Inf(-1)

	// compute scores and backtracking for remaining nodes
	var px int
	var sMax float64
	var bMax bi
	for x := 1; x < len(sg1); x++ {
		// compute sg1[x], "upper level"
		if x%stride == 0 {
			sg1[x] = math.Inf(-1) // sg1 -inf at "left edge"
		} else {
			// compute sg1[x], bg1[x] as max of two cases:
			px = x - g1
			// case 1: extend g1
			sMax = sg1[px] - gapExtendPenalty
			bMax = bi{g1, px}
			// case 2: open g1
			if s0 := smm[px] - gapOpenPenalty; s0 > sMax {
				sMax = s0
				bMax.level = mm
			}
			sg1[x] = sMax
			bg1[x] = bMax
		}

		// compute sg2[x], "lower level"
		if x < stride {
			sg2[x] = math.Inf(-1)
		} else {
			// compute sg2[x], bg2[x] as max of two cases:
			px = x - g2
			// case 1: extend g2
			sMax = sg2[px] - gapExtendPenalty
			bMax = bi{g2, px}
			// case 2: open g2
			if s0 := smm[px] - gapOpenPenalty; s0 > sMax {
				sMax = s0
				bMax.level = mm
			}
			sg2[x] = sMax
			bg2[x] = bMax
		}

		// compute smm[x], bmm[x] as max of three cases:
		// case 1:  close gap in s2
		sMax = sg2[x]
		bMax = bi{g2, x}
		// case 2: match/mismatch
		i := x / stride
		j := x % stride
		if i > 0 && j > 0 {
			px := x - mm
			fmt.Println("x, i, j:", x, i, j)
			s0 := smm[px] + float64(a.Score(s1[i-1], s2[j-1]))
			if s0 > sMax {
				sMax = s0
				bMax = bi{mm, px}
			}
		}
		// case 3: close gap in s1
		if s0 := sg1[x]; s0 > sMax {
			sMax = s0
			bMax = bi{g1, x}
		}
		// store accumulated max of three cases
		smm[x] = sMax
		bmm[x] = bMax
	}
	fmt.Println("final scores")
	for x := range sg1 {
		fmt.Printf("(%.0f %.0f %.0f) ", sg1[x], smm[x], sg2[x])
		if x%stride == stride-1 {
			fmt.Println()
		}
	}
	bl := map[int]string{g1: "g1", g2: "g2", mm: "mm"}
	fmt.Println("final bg2")
	for x, bx := range bg2 {
		px := bx.x
		fmt.Printf("(%s %d %d) ", bl[bx.level], px/stride, px%stride)
		if x%stride == stride-1 {
			fmt.Println()
		}
	}
	fmt.Println("final bmm")
	for x, bx := range bmm {
		px := bx.x
		fmt.Printf("(%s %d %d) ", bl[bx.level], px/stride, px%stride)
		if x%stride == stride-1 {
			fmt.Println()
		}
	}
	fmt.Println("final bg1")
	for x, bx := range bg1 {
		px := bx.x
		fmt.Printf("(%s %d %d) ", bl[bx.level], px/stride, px%stride)
		if x%stride == stride-1 {
			fmt.Println()
		}
	}
	fmt.Println("read out")
	// read out from last position
	x := len(sg1) - 1
	score = smm[x]
	fmt.Println("final score:", score)
	b := bmm
	// backtrack (from x at last position) then reverse
	for ; x > 0; x = b[x].x {
		switch b[x].level {
		case g1:
			t1 = append(t1, GapSymbol)
			t2 = append(t2, s2[x%stride-1])
			b = bg1
		case g2:
			t1 = append(t1, s1[x/stride-1])
			t2 = append(t2, GapSymbol)
			b = bg2
		default:
			t1 = append(t1, s1[x/stride-1])
			t2 = append(t2, s2[x%stride-1])
			b = bmm
		}
		fmt.Println("x, t1, t2, b[x].x:", x, t1, t2, b[x].x)
	}
	last := len(t1) - 1
	for i := range t1[:len(t1)/2] {
		t1[i], t1[last-i] = t1[last-i], t1[i]
		t2[i], t2[last-i] = t2[last-i], t2[i]
	}
	return
}

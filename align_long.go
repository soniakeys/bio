package bio

import "strings"

// AlignLong does fragment assembly of long reads.
//
// Greedy algorithm,
// http://en.wikipedia.org/wiki/Sequence_assembly#Greedy_algorithm.
//
// All fragments must be of the same length.  Overlaps of less than half
// the fragment length are not considered, thus, it is a requirement that
// complete assembly be possible considering only overlaps of at least
// half the fragment length.
//
// If all fragments can be assembled into a single string, that string is
// returned, otherwise the empty string is returned.
func AlignLong(ss []string) string {
	switch len(ss) {
	case 0:
		return ""
	case 1:
		return ss[0]
	}
	// from the problem: reads are same length and overlap by at least half.
	half := len(ss[0]) / 2

	// greedy algorithm,
	// http://en.wikipedia.org/wiki/Sequence_assembly#Greedy_algorithm.
	// 1. calculate pairwise alignments of all fragments,
	type overlap struct {
		si, ti int // indexes into ss
		x      int // x + half = overlap
	}
	var overlaps []overlap
	for si, s := range ss {
		sSuf := s[len(s)-half:]
		for ti, t := range ss {
			if t == s {
				continue
			}
			if x := strings.Index(t, sSuf); x >= 0 {
				overlaps = append(overlaps, overlap{si, ti, x})
			}
		}
	}
	var left int
	for nFragments := len(ss); nFragments > 1; nFragments-- {
		// 2. choose two fragments with the largest overlap.
		xMax := -1
		var oMax int // index into overlaps
		for i, o := range overlaps {
			s := ss[o.si]
			t := ss[o.ti]
			if s == "" || t == "" || t == s {
				continue
			}
			if o.x > xMax {
				xMax = o.x
				oMax = i
			}
		}
		if xMax == -1 {
			return ""
		}
		// 3. merge chosen fragments
		left = overlaps[oMax].si
		right := overlaps[oMax].ti
		ss[left] += ss[right][xMax+half:]
		ss[right] = ""
		// update overlap list.  overlaps from si can be deleted, overlaps
		// from ti need to be updated to be overlaps from si
		for i := 0; i < len(overlaps); {
			switch {
			case overlaps[i].si == left:
				// delete by swapping with last
				last := len(overlaps) - 1
				if i < last {
					overlaps[i] = overlaps[last]
				}
				overlaps = overlaps[:last]
				continue
			case overlaps[i].si == right:
				overlaps[i].si = left
			}
			i++
		}
		// 4. repeat step 2. and 3. until only one fragment is left
	}
	return ss[left]
}

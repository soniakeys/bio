// Copyright 2014 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import "github.com/soniakeys/multiset"

// SolvePartialDigest solves the "partial digest problem," useful for example
// in restriction mapping.
//
// For a list of n distances, argument l must be a full-information
// partial digest with cardinality n(n-1)/2.  The element type must be int.
//
// SolvePartialDigest returns a list of possible solutions for distance lists.
func SolvePartialDigest(l multiset.Multiset) (p [][]int) {
	maxL := func() (max int) {
		for x := range l {
			if x := x.(int); x > max {
				max = x
			}
		}
		return
	}
	width := maxL()
	l.AddElementCount(width, -1)
	x := []int{0, width}
	popX := func() int {
		last := len(x) - 1
		n := x[last]
		x = x[:last]
		return n
	}
	delta := func(y int) multiset.Multiset {
		m := multiset.Multiset{}
		for _, x1 := range x {
			if d := y - x1; d < 0 {
				m[-d]++
			} else {
				m[d]++
			}
		}
		return m
	}
	var place func()
	place1 := func(y int) {
		if d := delta(y); multiset.Subset(d, l) {
			x = append(x, y)
			l.SubtractCounts(d)
			place()
			l.AddCounts(d)
			popX()
		}
	}
	place = func() {
		if len(l) == 0 {
			p = append(p, append([]int{}, x...))
			return
		}
		y := maxL()
		place1(y)
		place1(width - y)
	}
	place()
	l.AddElementCount(width, 1)
	return
}

// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"fmt"
	"sort"
	"strconv"
	"strings"
)

// aa.go
//
// Amino acid definitions

// AA type holds amino acid sequences.
type AA []byte

// AA20 type holds amino acid sequences.  Content should be strictly limited
// to the twenty proteinogenic symbols, in upper case.
type AA20 []byte

// AA20Alphabet is the set of 20 proteinogenic amino acid symbols.
const AA20Alphabet = "ACDEFGHIKLMNPQRSTVWY"

func init() {
	initIntegerMass()
}

// String returns an AA converted to a string.
func (s AA) String() string {
	return string(s)
}

// String returns an AA20 converted to a string.
func (s AA20) String() string {
	return string(s)
}

// AAIntegerMassTable holds masses in an array that can be indexed
// efficiently.
//
// Normally you should use the function AAIntegerMass rather than
// access the table directly.
//
// Note that the table length is 25 and holds 20 amino acids, so there
// are a few holes.
var AAIntegerMassTable [25]int

func initIntegerMass() {
	AAIntegerMassTable['A'-'A'] = 71

	AAIntegerMassTable['C'-'A'] = 103
	AAIntegerMassTable['D'-'A'] = 115
	AAIntegerMassTable['E'-'A'] = 129
	AAIntegerMassTable['F'-'A'] = 147
	AAIntegerMassTable['G'-'A'] = 57
	AAIntegerMassTable['H'-'A'] = 137
	AAIntegerMassTable['I'-'A'] = 113

	AAIntegerMassTable['K'-'A'] = 128
	AAIntegerMassTable['L'-'A'] = 113
	AAIntegerMassTable['M'-'A'] = 131
	AAIntegerMassTable['N'-'A'] = 114

	AAIntegerMassTable['P'-'A'] = 97
	AAIntegerMassTable['Q'-'A'] = 128
	AAIntegerMassTable['R'-'A'] = 156
	AAIntegerMassTable['S'-'A'] = 87
	AAIntegerMassTable['T'-'A'] = 101

	AAIntegerMassTable['V'-'A'] = 99
	AAIntegerMassTable['W'-'A'] = 186

	AAIntegerMassTable['Y'-'A'] = 163
}

// AAIntegerMass returns integer mass for one of the 20 proteinogenic amino
// acids.
//
// It panics if aa is not in range 'A'-'Y'.
func AAIntegerMass(aa byte) int {
	return AAIntegerMassTable[aa-'A']
}

// NumSubPepLinear computes the theoretical number of subpeptides of a
// linear peptide of length l.
func NumSubPepLinear(l int) int {
	return (l*(l+1))/2 + 1
}

// NumSubPepCyclic computes the theoretical number of subpeptides of a
// cyclic peptide of length l.
func NumSubPepCyclic(l int) int {
	return l * (l - 1)
}

// IntSpec represents a spectrum of integer masses
type IntSpec []int

// Equal tests element-wise equality.
func (x IntSpec) Equal(y IntSpec) bool {
	if len(x) != len(y) {
		return false
	}
	for i, m := range x {
		if y[i] != m {
			return false
		}
	}
	return true
}

// Amino acid integer masses.  Like type AA, represents a peptide,
// but coded by integer mass rather than symbol.
type AAInt []byte

// String returns a string of masses separated with '-' characters.
func (p AAInt) String() string {
	s := make([]string, len(p))
	for i, m := range p {
		s[i] = strconv.Itoa(int(m))
	}
	return strings.Join(s, "-")
}

// CyclicSpec synthesizes the theoretical integer spectrum for a cyclic peptide.
func (p AAInt) CyclicSpec() IntSpec {
	sum := 0
	s := IntSpec{sum}
	for i := range p {
		sum = 0
		for j := (i + 1) % len(p); j != i; j = (j + 1) % len(p) {
			sum += int(p[j])
			s = append(s, sum)
		}
	}
	s = append(s, sum+int(p[len(p)-1]))
	sort.Ints(s)
	return s
}

// The 18 uinique integer masses out of the 20 proteinogenic amino acids.
var am18 = AAInt{57, 71, 87, 97, 99, 101, 103, 113, 114,
	115, 128, 129, 131, 137, 147, 156, 163, 186}

// expand grows each peptide in l by appending a single amino acid to the end.
// it does this for all 18 unique amino acids.  the length of the result list e
// will be 18 times the length of the argument list l.
func expand(l []AAInt) (e []AAInt) {
	for _, p := range l {
		for _, a := range am18 {
			e = append(e, append(p[:len(p):len(p)], a))
		}
	}
	return
}

// IntSpecCounts represents an integer mass spectrum as counts of each
// unique mass.
//
// This data structure is termed a multiset.
//
// The map key is integer mass, the map value is the number of occurrences.
type IntSpecCounts map[int]int

// Counts builds an IntSpecCounts object corresponding to the receiver.
func (s IntSpec) Counts() IntSpecCounts {
	sm := IntSpecCounts{}
	for _, m := range s {
		sm[m]++
	}
	return sm
}

// returns true if peptide p is consistent with spectrum spec.
// spec is represented as a multiset of masses in the spectrum.
func (p AAInt) consistent(spec IntSpecCounts) bool {
	ps := IntSpecCounts{}
	for i := range p {
		sum := 0
		for j := i; j < len(p); j++ {
			sum += int(p[j])
			f := ps[sum] + 1
			if f > spec[sum] {
				return false
			}
			ps[sum] = f
		}
	}
	return true
}

// SeqCyclicTheo sequences a cyclic peptide from its theoretical integer
// spectrum.  That is, a spectrum that is correct and complete.
func (s IntSpec) SeqCyclicTheo() (r []AAInt) {
	smap := s.Counts()
	l := []AAInt{{}}
	for len(l) > 0 {
		l = expand(l)
		for i := 0; i < len(l); {
			p := l[i] // peptide in list
			switch {
			case p.CyclicSpec().Equal(s):
				r = append(r, p)
			case !p.consistent(smap):
			default:
				i++
				continue
			}
			last := len(l) - 1
			l[i] = l[last]
			l = l[:last]
		}
	}
	return r
}

func (p AAInt) score(target IntSpecCounts) (score int) {
	for m, n := range p.CyclicSpec().Counts() {
		if t := target[m]; n < t {
			score += n
		} else {
			score += t
		}
	}
	return
}

func (pl *lbd) cut(target IntSpecCounts, n int) {
	l := *pl
	if len(l) <= n {
		return
	}
	for i, c := range l {
		l[i].s = c.score(target)
	}
	sort.Sort(l)
	ns := l[n].s
	i := n + 1
	for l[i].s == ns {
		i++
		if i == len(l) {
			return
		}
	}
	*pl = l[:i]
}

type cand struct {
	AAInt
	m int // total mass in AAInt
	s int // score
}

type lbd []cand

func (l lbd) Len() int           { return len(l) }
func (l lbd) Less(i, j int) bool { return l[i].s > l[j].s }
func (l lbd) Swap(i, j int)      { l[i], l[j] = l[j], l[i] }

// expand grows each peptide in l by appending a single amino acid to the end.
// it does this for all 18 unique amino acids.  the length of the result list e
// will be 18 times the length of the argument list l.
func (pl *lbd) expand() {
	l := *pl
	e := lbd{}
	for _, c := range l {
		n := len(c.AAInt)
		p := c.AAInt[:n:n]
		for _, a := range am18 {
			e = append(e, cand{
				AAInt: append(p, a),
				m:     c.m + int(a)})
		}
	}
	*pl = e
}

// SeqCyclicExp sequences a cyclic peptide from an experimental integer
// spectrum.  That is, a spectrum with errors and omissions.
func (s IntSpec) SeqCyclicExp(n int) (r []AAInt) {
	pm := 0
	smap := s.Counts()
	for m := range smap {
		if m > pm {
			pm = m
		}
	}
	l := lbd{{}}
	max := 0
	for len(l) > 0 {
		l.expand()
		for i := 0; i < len(l); {
			c := l[i] // candidate in list
			switch {
			case c.m < pm:
				i++
				continue
			case c.m == pm:
				switch s := c.score(smap); {
				case s > max:
					r = []AAInt{c.AAInt}
					max = s
					fmt.Print("\nmax:", max)
					fmt.Print(".")
				case s == max:
					r = append(r, c.AAInt)
					fmt.Print(".")
				}
			}
			last := len(l) - 1
			l[i] = l[last]
			l = l[:last]
		}
		l.cut(smap, n)
	}
	fmt.Println()
	return
}

// Convolve returns the convolution of an integer spectrum.  That is, a list
// of all pairwise differences.
func (s IntSpec) Convolve() []int {
	c := make([]int, (len(s)*(len(s)-1))/2)
	i := 0
	for j := 1; j < len(s); j++ {
		m1 := s[j]
		for _, m2 := range s[:j] {
			switch d := m2 - m1; {
			case d < 0:
				d = -d
				fallthrough
			case d > 0:
				c[i] = d
				i++
			}
		}
	}
	return c[:i]
}

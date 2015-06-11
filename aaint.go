// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"sort"
	"strconv"
	"strings"

	"github.com/soniakeys/multiset"
)

// aaint.go
//
// Functions on amino acid integer masses

// AA20IntegerMassTable holds masses in an array that can be indexed
// efficiently.
//
// Normally you should use the function AA20IntegerMass rather than
// access the table directly.
//
// Note that the table length is 25 and holds 20 amino acids, so there
// are a few holes with zero values.
var AA20IntegerMassTable [25]uint8

func init() {
	AA20IntegerMassTable['A'-'A'] = 71

	AA20IntegerMassTable['C'-'A'] = 103
	AA20IntegerMassTable['D'-'A'] = 115
	AA20IntegerMassTable['E'-'A'] = 129
	AA20IntegerMassTable['F'-'A'] = 147
	AA20IntegerMassTable['G'-'A'] = 57
	AA20IntegerMassTable['H'-'A'] = 137
	AA20IntegerMassTable['I'-'A'] = 113

	AA20IntegerMassTable['K'-'A'] = 128
	AA20IntegerMassTable['L'-'A'] = 113
	AA20IntegerMassTable['M'-'A'] = 131
	AA20IntegerMassTable['N'-'A'] = 114

	AA20IntegerMassTable['P'-'A'] = 97
	AA20IntegerMassTable['Q'-'A'] = 128
	AA20IntegerMassTable['R'-'A'] = 156
	AA20IntegerMassTable['S'-'A'] = 87
	AA20IntegerMassTable['T'-'A'] = 101

	AA20IntegerMassTable['V'-'A'] = 99
	AA20IntegerMassTable['W'-'A'] = 186

	AA20IntegerMassTable['Y'-'A'] = 163
}

// AA20IntegerMass returns integer mass for one of the 20 proteinogenic amino
// acids.
//
// It panics if aa is not in range 'A'-'Y'.
func AA20IntegerMass(aa byte) uint8 {
	return AA20IntegerMassTable[aa-'A']
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

// MassCounts represents mass counts of an IntSpec or convolution, for example.
//
// Elements of the embedded multiset have the dynamic type int.
type MassCounts struct {
	multiset.Multiset
}

// Int returns receiver p as a sequence of integer masses.
func (p AA20) Int() AAInt {
	m := make(AAInt, len(p))
	for i, a := range p {
		m[i] = AA20IntegerMass(a)
	}
	return m
}

// Amino acid integer masses.  Like type AA, represents a peptide,
// but coded by integer mass rather than symbol.
type AAInt []uint8

// String returns a string of masses separated with '-' characters.
func (p AAInt) String() string {
	s := make([]string, len(p))
	for i, m := range p {
		s[i] = strconv.Itoa(int(m))
	}
	return strings.Join(s, "-")
}

// PrefixSpec constructs a prefix spectrum for a linear peptide.
//
// The empty string is not counted.  The result spectrum is of the
// same length as the input peptide p.
func (p AAInt) PrefixSpec() IntSpec {
	s := make(IntSpec, len(p))
	sum := 0
	for i, m := range p {
		sum += int(m)
		s[i] = sum
	}
	return s
}

// IdealSpec constructs a prefix-suffix spectrum for a linear peptide.
//
// The empty string is counted once and the full string is counted once,
// giving a spectrum length of exactly 2* the length of the input peptide p.
func (p AAInt) IdealSpec() IntSpec {
	s := make(IntSpec, 2*len(p))
	sum := 0
	for i, m := range p {
		s[i] = sum
		sum += int(m)
	}
	for i, m := range p {
		s[len(p)+i] = sum
		sum -= int(m)
	}
	sort.Ints(s)
	return s
}

// LinearSpec synthesizes the theoretical integer spectrum for a linear peptide.
func (p AAInt) LinearSpec() IntSpec {
	s := make(IntSpec, NumSubPepLinear(len(p)))
	i := 1
	for j := len(p) - 1; j >= 0; j-- {
		sum := 0
		for _, m := range p[j:] {
			sum += int(m)
			s[i] = sum
			i++
		}
	}
	sort.Ints(s)
	return s
}

// CyclicSpec synthesizes the theoretical integer spectrum for a cyclic peptide.
func (p AAInt) CyclicSpec() IntSpec {
	// TODO improve efficiency
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

var (
	// The 18 uinique integer masses out of the 20 proteinogenic amino acids.
	AA18Int AAInt

	// look up IUPAC symbol for one of 18 integer masses.  value is '-' for
	// all values in in AA18Int.
	AAFrom18 [AA18IntMax + 1]byte
)

// Max integer mass is 186 (for Tryptophan.)
const AA18IntMax = 186

func init() {
	AA18Int = AAInt{57, 71, 87, 97, 99, 101, 103, 113, 114,
		115, 128, 129, 131, 137, 147, 156, 163, 186}

	// populate AAFrom18
	for m := range AAFrom18 {
		AAFrom18[m] = '-'
	}
	for z := len(AA20IntegerMassTable) - 1; z >= 0; z-- {
		if m := AA20IntegerMassTable[z]; m > 0 {
			AAFrom18[m] = 'A' + byte(z)
		}
	}
}

func (p AAInt) AA20() (p2 AA20, ok bool) {
	p2 = make(AA20, len(p))
	ok = true
	for i, m := range p {
		s := AAFrom18[m]
		if s == '-' {
			ok = false
		}
		p2[i] = s
	}
	return
}

// expand grows each peptide in l by appending a single amino acid to the end.
// it does this for all 18 unique amino acids.  the length of the result list e
// will be 18 times the length of the argument list l.
func expand(l []AAInt) (e []AAInt) {
	// TODO allocate 18x up front
	for _, p := range l {
		for _, a := range AA18Int {
			e = append(e, append(p[:len(p):len(p)], a))
		}
	}
	return
}

// Counts constructs an MassCounts object representing the counts of
// distinct masses in the receiver IntSpec s.
func (s IntSpec) Counts() MassCounts {
	sm := multiset.Multiset{}
	for _, m := range s {
		sm[m]++
	}
	return MassCounts{sm}
}

// returns true if peptide p is consistent with spectrum spec.
// spec is represented as its counts.
func (p AAInt) consistent(spec MassCounts) bool {
	// TODO just a map would be fine here.
	ps := multiset.Multiset{}
	for i := range p {
		sum := 0
		for j := i; j < len(p); j++ {
			sum += int(p[j])
			f := ps[sum] + 1
			if f > spec.Multiset[sum] {
				return false
			}
			ps[sum] = f
		}
	}
	return true
}

// SeqCyclicTheo sequences a cyclic peptide from its theoretical integer
// spectrum.  That is, a spectrum that is correct and complete.
//
// Peptides are returned as integer mass sequences, drawn from the set of
// 18 unique masses of the 20 proteinogenic amino acids.  Returned is a list
// of all peptides with theoretical spectra consistent with receiver s.
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

// CyclicCommonCounts returns the number of mass counts in common between a
// cyclic peptide pep and a spectrum.
func (pep AAInt) CyclicCommonCounts(spec MassCounts) int {
	return commonCounts(pep.CyclicSpec(), spec)
}

// LinearCommonCounts returns the number of mass counts in common between a
// linear peptide and a spectrum.
func (pep AAInt) LinearCommonCounts(spec MassCounts) int {
	return commonCounts(pep.LinearSpec(), spec)
}

func commonCounts(spec IntSpec, specCounts MassCounts) int {
	return multiset.IntersectionCardinality(spec.Counts().Multiset, specCounts.Multiset)
}

// SeqCyclic18Exp sequences a cyclic peptide from an experimental integer
// spectrum.  That is, a spectrum with errors and omissions.
//
// The amino acid integer mass set used is the 18 unique masses of the 20
// proteinogenic amino acids.
//
// Argument n is a kind of a search width.  Small numbers may miss solutions.
// Large numbers waste time.  Try a number on the order of the parent mass.
//
// Peptides are returned as integer mass sequences.
// Returned peptides will have mass matching the parent mass of s, and be those
// with spectra best matching s.  All peptides "tied" for the best match are
// returned.
func (s IntSpec) SeqCyclic18Exp(n int) []AAInt {
	return s.leaderboard(n, 1, AA18Int)
}

/*
func (s IntSpec) SeqCyclic200Exp(n int) []AAInt {
	var a AAInt
	for i := uint8(57); i <= 200; i++ {
		a = append(a, i)
	}
	return s.leaderboard(n, 1, a)
}
*/

// SeqCyclicExp sequences a cyclic peptide from an experimental integer
// spectrum.  That is, a spectrum with errors and omissions.
//
// The amino acid integer masses are not restricted to those of the 20
// proteinogenic amino acids.  The mass set is picked with a convolution
// technique.
//
// Argument n is a kind of a search width.  Small numbers may miss solutions.
// Large numbers waste time.  Try a number on the order of the parent mass.
//
// Peptides are returned as integer mass sequences.
// Returned peptides will have mass matching the parent mass of s, and be those
// with spectra best matching s.
//
// Argument nr is a target number of results to return.  If more than nr
// peptides are found matching the parent mass, nr are returned plus any
// that are "tying" for last place.
func (s IntSpec) SeqCyclicExp(alphabetLen, n, nr int) (r []AAInt) {
	c := s.Convolve()
	a := c.cutAA(alphabetLen)
	return s.leaderboard(n, nr, a)
}

func (s IntSpec) leaderboard(n, nr int, set AAInt) []AAInt {
	pm := 0 // parent mass (code here doesn't count on it being s[len(s)-1]
	smap := s.Counts()
	for m := range smap.Multiset {
		if m := m.(int); m > pm {
			pm = m
		}
	}
	l := lbd{{}}  // leaderboard, initialized with 0-peptide
	done := lbd{} // candidates matching parent mass
	for len(l) > 0 {
		l.expand(set)
		for i := 0; i < len(l); {
			c := l[i] // candidate in list
			switch {
			case c.m == pm:
				c.s = c.AAInt.CyclicCommonCounts(smap) // final score
				done = append(done, c)
				fallthrough
			case c.m < pm:
				i++
				continue // candidate remains in the running
			}
			// c.m >= pm: remove from leaderboard
			last := len(l) - 1
			l[i] = l[last]
			l = l[:last]
		}
		if n < len(l) {
			for i, c := range l {
				l[i].s = c.AAInt.LinearCommonCounts(smap)
			}
			l = l[:Cut(l, n)]
		}
	}
	nr = Cut(done, nr)
	r := make([]AAInt, nr)
	for i, c := range done[:nr] {
		r[i] = c.AAInt
	}
	return r
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

// expand grows each peptide in receiver leaderboard pl by appending a single
// amino acid mass to the end.  it does this for each amino acid mass in
// argument set.  the length of the resulting leader board will be the
// original length times the length of the argument list l.
func (pl *lbd) expand(set AAInt) {
	l := *pl   // original leaderboard
	e := lbd{} // new leaderboard to be populated
	for _, c := range l {
		n := len(c.AAInt)
		p := c.AAInt[:n:n] // copy the peptide
		for _, a := range set {
			e = append(e, cand{
				AAInt: append(p, a),
				m:     c.m + int(a)})
		}
	}
	*pl = e // replace original leaderboard with new expanded leaderboard
}

// Cut finds the first n values in list, including ties for last place.
//
// The list are argument is sorted as a side effect of this function.
// There is no need for it to already be sorted.
//
// Returned is the index for slicing the now-sorted list to cut at the first
// n values.
func Cut(list sort.Interface, n int) int {
	if n <= 0 {
		return 0
	}
	ln := list.Len()
	sort.Sort(list)
	if ln <= n {
		return ln
	}
	r := n
	for nth := n - 1; r < ln && !list.Less(nth, r); r++ {
	}
	return r
}

// Convolve returns the convolution of an integer spectrum.
//
// The convolution consists of all positive pairwise differences between
// elements of receiver IntSpec s.  Distinct mass differences are counted
// and returned as the result.
func (s IntSpec) Convolve() MassCounts {
	m := multiset.Multiset{}
	for _, m1 := range s {
		for _, m2 := range s {
			if d := m1 - m2; d > 0 {
				m[d]++
			}
		}
	}
	return MassCounts{m}
}

func (m MassCounts) cutAA(cut int) AAInt {
	var cand mfc
	for mass, freq := range m.Multiset {
		if m := mass.(int); m >= 57 && m <= 200 {
			cand = append(cand, mf{m, freq})
		}
	}
	cut = Cut(cand, cut)
	r := make(AAInt, cut)
	for i, mf := range cand[:cut] {
		r[i] = uint8(mf.mass)
	}
	return r
}

type mf struct {
	mass int
	freq int
}

type mfc []mf

func (m mfc) Len() int           { return len(m) }
func (m mfc) Less(i, j int) bool { return m[i].freq > m[j].freq }
func (m mfc) Swap(i, j int)      { m[i], m[j] = m[j], m[i] }

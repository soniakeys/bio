// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"fmt"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/soniakeys/multiset"
)

// aa.go
//
// Amino acid definitions

// AA type holds amino acid sequences.
//
// Symbols are not restricted to AA20Alphabet, data may include other
// symbols such a gap symbols.
type AA []byte

// String returns an AA converted to a string.
func (s AA) String() string {
	return string(s)
}

// AA20 type holds amino acid sequences.  Content should be strictly limited
// to the symbols of AA20Alphabet, in upper case.  Methods may panic
// on other symbols.
type AA20 []byte

var AA20Alphabet = AA20("ACDEFGHIKLMNPQRSTVWY") // IUPAC symbols for the 20 proteinogenic amino acids

// String returns an AA20 converted to a string.
func (s AA20) String() string {
	return string(s)
}

// AA20IntegerMassTable holds masses in an array that can be indexed
// efficiently.
//
// Normally you should use the function AA20IntegerMass rather than
// access the table directly.
//
// Note that the table length is 25 and holds 20 amino acids, so there
// are a few holes.
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

// MassCounts represents mass counts of an IntSpec or convolution, for example.
//
// Elements of the embedded multiset have the dynamic type int.
type MassCounts struct {
	multiset.Multiset
}

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
var AA18Int = AAInt{57, 71, 87, 97, 99, 101, 103, 113, 114,
	115, 128, 129, 131, 137, 147, 156, 163, 186}

// expand grows each peptide in l by appending a single amino acid to the end.
// it does this for all 18 unique amino acids.  the length of the result list e
// will be 18 times the length of the argument list l.
func expand(l []AAInt) (e []AAInt) {
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

func (pep AAInt) score(spec MassCounts) int {
	pc := pep.CyclicSpec().Counts()
	return multiset.IntersectionCardinality(pc.Multiset, spec.Multiset)
}

func (pl *lbd) cut(target MassCounts, n int) {
	l := *pl
	if len(l) <= n {
		return
	}
	for i, c := range l {
		l[i].s = c.AAInt.score(target)
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

// SeqCyclicExp sequences a cyclic peptide from an experimental integer
// spectrum.  That is, a spectrum with errors and omissions.
//
// The amino acid integer mass set used is the 18 unique masses of the 20
// proteinogenic amino acids.
//
// Peptides are returned as integet mass sequences.  Returned is a list of
// all peptides found with
func (s IntSpec) SeqCyclic18Exp(n int) []AAInt {
	return s.leaderboard(n, AA18Int)
}

func (s IntSpec) SeqCyclicExp(alphabetLen, n int) (r []AAInt) {
	c := s.Convolve()
	fmt.Println("c:", c)
	a := c.cutAA(alphabetLen)
	fmt.Println("a:", a)
	return s.leaderboard(n, a)
}

func (s IntSpec) leaderboard(n int, set AAInt) (r []AAInt) {
	pm := 0
	smap := s.Counts()
	for m := range smap.Multiset {
		if m := m.(int); m > pm {
			pm = m
		}
	}
	l := lbd{{}}
	max := 0
	for len(l) > 0 {
		l.expand(set)
		for i := 0; i < len(l); {
			c := l[i] // candidate in list
			switch {
			case c.m < pm:
				i++
				continue
			case c.m == pm:
				switch s := c.AAInt.score(smap); {
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

// AA20MonoisotopicMassTable holds masses in an array that can be indexed
// efficiently.
//
// Normally you should use the function AA20MonoisotopicMass rather than
// access the table directly.
//
// Note that the table length is 25 and holds 20 amino acids, so there
// are a few holes.
var AA20MonoisotopicMassTable [25]float64

// Note: must precede init of aaByMass.
func init() {
	AA20MonoisotopicMassTable['A'-'A'] = 71.03711

	AA20MonoisotopicMassTable['C'-'A'] = 103.00919
	AA20MonoisotopicMassTable['D'-'A'] = 115.02694
	AA20MonoisotopicMassTable['E'-'A'] = 129.04259
	AA20MonoisotopicMassTable['F'-'A'] = 147.06841
	AA20MonoisotopicMassTable['G'-'A'] = 57.02146
	AA20MonoisotopicMassTable['H'-'A'] = 137.05891
	AA20MonoisotopicMassTable['I'-'A'] = 113.08406

	AA20MonoisotopicMassTable['K'-'A'] = 128.09496
	AA20MonoisotopicMassTable['L'-'A'] = 113.08406
	AA20MonoisotopicMassTable['M'-'A'] = 131.04049
	AA20MonoisotopicMassTable['N'-'A'] = 114.04293

	AA20MonoisotopicMassTable['P'-'A'] = 97.05276
	AA20MonoisotopicMassTable['Q'-'A'] = 128.05858
	AA20MonoisotopicMassTable['R'-'A'] = 156.10111
	AA20MonoisotopicMassTable['S'-'A'] = 87.03203
	AA20MonoisotopicMassTable['T'-'A'] = 101.04768

	AA20MonoisotopicMassTable['V'-'A'] = 99.06841
	AA20MonoisotopicMassTable['W'-'A'] = 186.07931

	AA20MonoisotopicMassTable['Y'-'A'] = 163.06333
}

const (
	AA20Lightest = 'G' // Glycine
	AA20Heaviest = 'W' // Tryptophan
)

// aaMass associates an amino acid symbol with a mass in a struct.
type aaMass struct {
	aa   byte
	mass float64
}

// aaMassList is a container type satisfying sort.Interface.
type aaMassList []aaMass

func (a aaMassList) Len() int           { return len(a) }
func (a aaMassList) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a aaMassList) Less(i, j int) bool { return a[i].mass < a[j].mass }

// aaByMass holds a list of the 20 amino acids sorted by mass,
// in increasing order.
var aaByMass aaMassList

// Note:  this must follow init of AA20MonoisotopicMassTable.
func init() {
	aaByMass = make(aaMassList, len(AA20Alphabet))
	for i := 0; i < len(AA20Alphabet); i++ {
		aa := AA20Alphabet[i]
		aaByMass[i] = aaMass{aa, AA20MonoisotopicMass(aa)}
	}
	sort.Sort(aaByMass)
}

// AA20MonoisotopicMass returns the residue mass of an amino acid.
//
// The residue is the amino acid without its water group.
// Monoisotopic means it is the mass assuming the most common isotope
// of each element.  The unit is daltons.
//
// Symbols not in AA20Alphabet may cause panic.
func AA20MonoisotopicMass(aa byte) float64 {
	return AA20MonoisotopicMassTable[aa-'A']
}

// AA20NearestMass takes a mass m and returns the symbol of the amino acid with
// monoisotopic residue mass nearest m.  For convenience, it also returns the
// absolute value of the mass difference.
func AA20NearestMass(m float64) (aa byte, d float64) {
	i := sort.Search(20, func(i int) bool { return aaByMass[i].mass >= m })
	switch i {
	case 0:
	case 20:
		i--
	default:
		di := aaByMass[i].mass - m
		dLess := m - aaByMass[i-1].mass
		if dLess < di {
			i--
		}
	}
	aam := aaByMass[i]
	return aam.aa, math.Abs(m - aam.mass)
}

const (
	WaterMassAverage      = 18.01528 // sum of constituent average atomic masses
	WaterMassMonoisotopic = 18.01056 // sum of constituent monoisotopic masses
)

// Weight returns the monoisotopic weight of the amino acid chain.
func (a AA20) Weight() float64 {
	return a.sum(true)
}

// perhaps overkill, but the divide and conquer algorithm should preserve
// precision a little better than linear summing.
func (a AA20) sum(addWater bool) (s float64) {
	if len(a) < 10 {
		if addWater {
			s = WaterMassMonoisotopic
		}
		for _, aa := range a {
			s += AA20MonoisotopicMass(aa)
		}
		return s
	}
	half := len(a) / 2
	return a[:half].sum(addWater) + a[half:].sum(false)
}

// PrefixSpectrum returns a monoisotopic spectrum of all prefixes of the
// amino acid string.
func (s AA20) PrefixSpectrum() []float64 {
	r := make([]float64, len(s)+1)
	for i, aa := range s {
		r[i+1] = r[i] + AA20MonoisotopicMass(aa)
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
	if len(cand) < cut {
		cut = len(cand)
	} else {
		sort.Sort(cand)
		fcut := cand[cut-1].freq
		for cut < len(cand) && cand[cut].freq == fcut {
			cut++
		}
	}
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

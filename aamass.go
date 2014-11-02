// Copyright 2013 Sonia Keys
// License MIT: http://opensource.org/licenses/MIT

package bio

import (
	"fmt"
	"math"
	"sort"
)

// aamass.go
//
// Functions on amino acid masses as float64s.
// (Integer mass functions are in aaint.go.)

// AA20MonoisotopicMassTable holds masses in an array that can be indexed
// efficiently.
//
// Normally you should use the function AA20MonoisotopicMass rather than
// access the table directly.
//
// Note that the table length is 25 and holds 20 amino acids, so there
// are a few holes.
var AA20MonoisotopicMassTable [25]float64

// Note: must precede init of aa20ByMass.
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

const ProtonMass = 1.007

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

// aa20ByMass holds a list of the 20 amino acids sorted by mass,
// in increasing order.
var aa20ByMass aaMassList

// Note:  this must follow init of AA20MonoisotopicMassTable.
func init() {
	aa20ByMass = make(aaMassList, len(AA20Alphabet))
	for i := 0; i < len(AA20Alphabet); i++ {
		aa := AA20Alphabet[i]
		aa20ByMass[i] = aaMass{aa, AA20MonoisotopicMass(aa)}
	}
	sort.Sort(aa20ByMass)
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
	i := sort.Search(20, func(i int) bool { return aa20ByMass[i].mass >= m })
	switch i {
	case 0:
	case 20:
		i--
	default:
		di := aa20ByMass[i].mass - m
		dLess := m - aa20ByMass[i-1].mass
		if dLess < di {
			i--
		}
	}
	aam := aa20ByMass[i]
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

// AAMass represents an amino acid sequence as a sequence of the masses of
// the individual amino acids.
type AAMass []float64

// MonoisotopicMass returns receiver p as a sequences of masses.
func (p AA20) MonoisotopicMass() AAMass {
	m := make(AAMass, len(p))
	for i, a := range p {
		m[i] = AA20MonoisotopicMass(a)
	}
	return m
}

// LinearSpec synthesizes the theoretical monoisotopic spectrum for
// a linear peptide.
func (p AAMass) LinearSpec() MassSpec {
	s := make(MassSpec, NumSubPepLinear(len(p)))
	i := 1
	for j := len(p) - 1; j >= 0; j-- {
		sum := 0.
		for _, m := range p[j:] {
			sum += m
			s[i] = sum
			i++
		}
	}
	sort.Float64s(s)
	return s
}

// CyclicSpec synthesizes the theoretical monoisotopic spectrum for
// a cyclic peptide.
func (p AAMass) CyclicSpec() MassSpec {
	// TODO improve efficiency
	sum := 0.
	s := MassSpec{sum}
	for i := range p {
		sum = 0
		for j := (i + 1) % len(p); j != i; j = (j + 1) % len(p) {
			sum += p[j]
			s = append(s, sum)
		}
	}
	s = append(s, sum+p[len(p)-1])
	sort.Float64s(s)
	return s
}

// PrefixSpectrum returns a monoisotopic spectrum of all prefixes of the
// amino acid string.
func (s AA20) PrefixSpectrum() MassSpec {
	r := make(MassSpec, len(s)+1)
	for i, aa := range s {
		r[i+1] = r[i] + AA20MonoisotopicMass(aa)
	}
	return r
}

// MassSpec represents a mass spectrum.  Generally the values are sorted
// and the last value represents a parent mass.
type MassSpec []float64

// NearestFloat64 finds the element of s nearest x.
//
// Slice s must be sorted in increasing order.
//
// Retuned is the index into s of the nearest value and the absolute value
// of the difference between x and the nearest value.
func NearestFloat64(s []float64, x float64) (i int, d float64) {
	i = sort.Search(len(s), func(i int) bool { return s[i] >= x })
	switch i {
	case 0:
	case len(s):
		i--
	default:
		di := s[i] - x
		dLess := x - s[i-1]
		if dLess < di {
			i--
		}
	}
	return i, math.Abs(x - s[i])
}

// score is number of masses in theo within .3 of a mass in measured
func (meas MassSpec) score(theo MassSpec) int {
	fmt.Println("theo:", theo)
	c := 0
	lastTm := 0.
	lastMatch := false
	for _, tm := range theo {
		if tm == lastTm {
			if lastMatch {
				c++
			}
			continue
		}
		_, d := NearestFloat64(meas, tm)
		fmt.Print(tm, " within ", d)
		lastTm = tm
		lastMatch = d <= .3
		if lastMatch {
			fmt.Println(" -- counted")
			c++
		} else {
			fmt.Println()
		}
	}
	fmt.Println("Score:", c)
	return c
}

// SeqCyclic20 sequences a cyclic peptide from an experimental spectrum.
//
// The amino acid mass set used is the monoisotopic masses of the 20
// proteinogenic amino acids.
//
// Argument n is a kind of a search width.  Small numbers may miss solutions.
// Large numbers waste time.  Try a number on the order of the parent mass.
//
// Returned peptides will have mass matching the parent mass of s, and be those
// with spectra best matching s.  All peptides "tied" for the best match are
// returned.
func (s MassSpec) SeqCyclic20(n int) []AA20 {
	l := s.leaderboard(n, 1)
	r := make([]AA20, len(l))
	for i, c := range l {
		r[i] = c.AA20
	}
	return r
}

func (s MassSpec) leaderboard(n, nr int) []massCand {
	sort.Float64s(s)
	pm := s[len(s)-1] // parent mass
	// leaderboard, initialized with 0-peptide
	l := massLbd{{}}
	done := massLbd{} // candidates matching parent mass
	for len(l) > 0 {
		l.expand20()
		for i := 0; i < len(l); {
			c := l[i] // candidate in list
			switch {
			case c.Mass < pm:
				i++
				continue // candidate remains in the running
			case math.Abs(pm-c.Mass) <= .3:
				c.Score = s.score(c.AAMass.CyclicSpec()) // final Score
				done = append(done, c)
			}
			// c.m >= pm: remove from leaderboard
			last := len(l) - 1
			l[i] = l[last]
			l = l[:last]
		}
		if n < len(l) {
			for i, c := range l {
				l[i].Score = s.score(c.AAMass.LinearSpec())
			}
			l = l[:Cut(l, n)]
		}
	}
	nr = Cut(done, nr)
	return done[:nr]
}

type massCand struct {
	AA20
	AAMass
	Mass  float64 // total mass in AA20
	Score int     // Score
}

type massLbd []massCand

func (l massLbd) Len() int           { return len(l) }
func (l massLbd) Less(i, j int) bool { return l[i].Score > l[j].Score }
func (l massLbd) Swap(i, j int)      { l[i], l[j] = l[j], l[i] }

func (pl *massLbd) expand20() {
	l := *pl       // original leaderboard
	e := massLbd{} // new leaderboard to be populated
	for _, c := range l {
		n := len(c.AA20)
		pep := c.AA20[:n:n] // force append to copy the peptide
		pepM := c.AAMass[:n:n]
		for _, a := range AA20Alphabet {
			m := AA20MonoisotopicMass(a)
			e = append(e, massCand{
				AA20:   append(pep, a),
				AAMass: append(pepM, m),
				Mass:   c.Mass + m})
		}
	}
	*pl = e // replace original leaderboard with new expanded leaderboard
}

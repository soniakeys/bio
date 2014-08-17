package bio

import (
	"math"
	"strconv"
)

// ParseFloats parses a string slice to a slice of float64s.
//
// Each string must contain a single slice trimmed of white space.
// If an error is encountered at any point, nil and the error are returned.
func ParseFloats(ss []string) (s []float64, err error) {
	s = make([]float64, len(ss))
	for i, m := range ss {
		if s[i], err = strconv.ParseFloat(m, 64); err != nil {
			return nil, err
		}
	}
	return
}

// Quantize converts a float64 slice to an int64 slice.
//
// Values are first multiplied by the scale parameter, then rounded to int64s.
func Quantize(s []float64, scale float64) []int64 {
	r := make([]int64, len(s))
	for i, m := range s {
		r[i] = int64(math.Floor(m*scale + .5))
	}
	return r
}

// ConvolveSectra returns a mass difference corresponding to the spectral
// convolution of two quantized spectra.
//
// The shift is computed with a Minkowski difference that produces a multiset
// of mass differences (mass2 - mass1) and the multiplicity of peaks of those
// mass differences.
// The mass difference is returned as diff, and the multiplicity as mult.
func ConvolveSpectra(spec1, spec2 []int64) (diff int64, mult int) {
	// ms is the spectral multiset representing the Minkowski difference
	ms := make(map[int64]int)
	// Minkowski difference algorithm
	for _, m1 := range spec1 {
		for _, m2 := range spec2 {
			d := m2 - m1
			ms[d]++
			n := ms[d]
			// Accumulate max multiplicity in the same pass.
			// This avoids a separate loop over the entire multiset after
			// construction.
			if n > mult {
				mult = n
				diff = d
			}
		}
	}
	return
}

package bio

// Mendel1 computes offspring genotype probabilities by Mendel's first
// law, the law of segregation.
//
// Inputs are the number individuals in a population which are dominant,
// heterozygous, and recessive respectively for a specific trait.
//
// Returned are the probabilities of an offspring being dominant,
// heterozygous, or recessive for the trait.
func Mendel1(dom, het, rec int) (pd, ph, pr float64) {
	dd := float64(dom * (dom - 1))
	hh := float64(het * (het - 1))
	rr := float64(rec * (rec - 1))
	dh := float64(dom * het)
	hr := float64(het * rec)
	rd := float64(rec * dom)
	pop := dom + het + rec // total population
	denom := float64(pop * (pop - 1))
	return (dd + dh + hh*.25) / denom,
		(rd*2 + dh + hr + hh*.5) / denom,
		(rr + hr + hh*.25) / denom
}

// Permute computes permutations of an integer slice.
//
// It takes a callback function that Permute calls on each permutation.
// The callback function should normally return false.  If it returns true,
// Permute terminates early.
func Permute(s []int, emit func(perm []int) (terminate bool)) {
	if len(s) == 0 {
		emit(s)
		return
	}
	// A fairly simple recursive double swap algorithm.
	// Not lexicographic order, not the fastest, but still pretty fast.
	a := append([]int{}, s...)
	var rc func(int) bool
	rc = func(last int) bool {
		if last == 0 {
			return emit(a)
		}
		for i := 0; i <= last; i++ {
			a[i], a[last] = a[last], a[i]
			if rc(last - 1) {
				return true
			}
			a[i], a[last] = a[last], a[i]
		}
		return false
	}
	rc(len(s) - 1)
}

// PermuteAll generates permutations of an integer slice.
//
// It is like Permute but without the early termination feature.
func PermuteAll(s []int, emit func(perm []int)) {
	if len(s) == 0 {
		emit(s)
		return
	}
	a := append([]int{}, s...)
	var rc func(int)
	rc = func(last int) {
		if last == 0 {
			emit(a)
		}
		for i := 0; i <= last; i++ {
			a[i], a[last] = a[last], a[i]
			rc(last - 1)
			a[i], a[last] = a[last], a[i]
		}
	}
	rc(len(s) - 1)
}

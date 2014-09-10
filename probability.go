package bio

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

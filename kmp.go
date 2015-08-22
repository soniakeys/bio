package bio

// KMP type represents a prepocessed pattern as used by the
// Knuth-Morris-Pratt string search algorithm.
//
// Member P is the Seq passed to the constructor.
// KMP member functions will not work if the Seq is subsequently modified.
// Copy the Seq if needed.
type KMP struct {
	P Seq
	b []int
}

// NewKMP creates and generates a KMP object from pattern p.
//
// Argument P is stored as a struct member.
// KMP member functions will not work if the Seq is subsequently modified.
// Copy the Seq if needed.
func NewKMP(P Seq) KMP {
	// ref: http://www.inf.fh-flensburg.de/lang/algorithmen/pattern/kmpen.htm
	b := make([]int, len(P)+1)
	i, j := 0, -1
	b[0] = j
	for i < len(P) {
		for j >= 0 && P[i] != P[j] {
			j = b[j]
		}
		i++
		j++
		b[i] = j
	}
	return KMP{P, b}
}

// Index returns the slice index of the first match of k.P in t.
// It returns -1 if there is no match.
func (k KMP) Index(t Seq) int {
	i, j := 0, -1
	for i < len(t) {
		for j >= 0 && t[i] != k.P[j] {
			j = k.b[j]
		}
		i++
		j++
		if j == len(k.P) {
			return i - j
		}
	}
	return -1
}

// AllIndex returns indexes of all matches of k.P in t.
func (k KMP) AllIndex(t Seq) (x []int) {
	i, j := 0, -1
	for i < len(t) {
		for j >= 0 && t[i] != k.P[j] {
			j = k.b[j]
		}
		i++
		j++
		if j == len(k.P) {
			x = append(x, i-j)
			j = k.b[j]
		}
	}
	return
}

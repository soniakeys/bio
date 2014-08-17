package bio

// KMP type represents a "failure array" as used by the Knuth-Morris-Pratt
// algorithm.
type KMP []int

// NewKMP creates and generates a KMP object from byte slice s.
func NewKMP(s []byte) KMP {
	k := make(KMP, len(s))
	if len(s) < 2 {
		return k
	}
	for pos, cnd := 1, 0; pos < len(s); {
		switch {
		case s[pos] == s[cnd]:
			cnd++
			k[pos] = cnd
			pos++
		case cnd > 0:
			cnd = k[cnd-1]
		default:
			pos++
		}
	}
	return k
}

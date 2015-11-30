package bio

import "sort"

// BWT is a searchable Burrows-Wheeler transform of a string.
type BWT struct {
	bwt []byte
	psa map[int]int // partial suffix array
	fo  [256]int    // first occurrence
	cx  [256]int    // index from byte to counts
	cp  [][]int     // checkpoints
	mod int         // sparseness of checkpoint array
}

// NewBWT constructs a Burrows-Wheeler transform of a string.
//
// The string is interpreted byte-wise, not as runes.
// The sentinal byte must be a byte with a value less than that of any byte
// present in the string.
// It can be present as the last byte of the string or can be missing.
// If missing, it will be added to an internal copy of the string (increasing
// memory requirements.)  The last parameter, mod, is a tunable parameter for
// reducing memory requirements.  For long strings try a value of 100.
func NewBWT(s string, sentinal byte, mod int) *BWT {
	if s[len(s)-1] != sentinal {
		s += string(sentinal)
	}
	bx := &BWT{
		bwt: make([]byte, len(s)),
		psa: map[int]int{},
		cp:  make([][]int, len(s)/mod+1),
		mod: mod,
	}
	// construct bwt, psa
	ss := make([]string, len(s))
	for i := range s {
		ss[i] = s[i:]
	}
	sort.Strings(ss)
	for i, si := range ss {
		sa := len(s) - len(si) // suffix array value
		if sa%mod == 0 {
			bx.psa[i] = sa // only store partial suffix array
		}
		bx.bwt[i] = s[(sa+len(s)-1)%len(s)]
	}

	lex := append(bwb{}, bx.bwt...)
	sort.Sort(lex)
	sx := []byte{} // unique symbols
	for i := range lex {
		b := lex[i]
		if bx.cx[b] == 0 {
			bx.fo[b] = i
			bx.cx[b] = len(sx)
			sx = append(sx, b)
		}
	}
	counts := make([]int, len(sx))
	for i, b := range bx.bwt {
		if i%mod == 0 {
			bx.cp[i/mod] = append([]int{}, counts...)
		}
		counts[bx.cx[b]]++
	}
	if len(bx.bwt)%mod == 0 {
		bx.cp[len(bx.cp)-1] = counts
	}
	return bx
}

func (b *BWT) count(sym byte, i int) int {
	m := i / b.mod
	c := b.cp[m][b.cx[sym]]
	for j := m * b.mod; j < i; j++ {
		if b.bwt[j] == sym {
			c++
		}
	}
	return c
}

// AllIndex returns all indexes of pat in the string indexed by BWT.
func (b *BWT) AllIndex(pat string) []int {
	top := 0
	bot := len(b.bwt) - 1
	for last := len(pat) - 1; top <= bot; last-- {
		if last < 0 {
			p := make([]int, bot+1-top)
			for i := range p {
				j := top + i
				for k := 0; ; k++ {
					m, ok := b.psa[j]
					if ok {
						p[i] = m + k
						break
					}
					sym := b.bwt[j]
					j = b.fo[sym] + b.count(sym, j)
				}
			}
			return p
		}
		sym := pat[last]
		top = b.fo[sym] + b.count(sym, top)
		bot = b.fo[sym] + b.count(sym, bot+1) - 1
	}
	return nil
}

type bwb []byte

func (l bwb) Len() int           { return len(l) }
func (l bwb) Swap(i, j int)      { l[i], l[j] = l[j], l[i] }
func (l bwb) Less(i, j int) bool { return l[i] < l[j] }

package bio

import (
	"errors"
	"math/big"
)

// CharacterArray represents a subset of phylogenic taxa that are present
// in some situation.
//
// The string representation is a Go string of '1's and '0's where a 1
// represents a taxon present and a 0 represents a taxon not present.
// The internal representation uses efficient bit-wise operations.
// Methods such as IncludeTaxon take indexes that correspond to indexes
// of the string representation.
type CharacterArray struct {
	n        int
	universe *big.Int
	bits     *big.Int
}

// NewCharacterArray constructs a CharacterArray object.
//
// Argument n is the number of taxa represented.
//
// Argument u is stored internally as a "universe" object.  The universe
// object is a bit string of length n of all 1's.
// Universe objects are constant after construction and so can be shared
// between CharacterArray objects.  A u argument that has a bit length of n
// will be reused and the allocation will be avoided.
// Otherwise a new universe object will be constructed.
//
// Argument s is an initial string.  See CharacterArray.SetString().
func NewCharacterArray(n int, u *big.Int, s string) *CharacterArray {
	ca := &CharacterArray{
		n:        n,
		universe: u,
		bits:     new(big.Int),
	}
	setU(ca)
	return ca.SetString(s)
}

var one = big.NewInt(1)

func setU(a *CharacterArray) {
	u := a.universe
	if u == nil || u.BitLen() != a.n {
		u = new(big.Int)
	}
	u.Sub(u.Lsh(one, uint(a.n)), one)
	a.universe = u
}

func setN(a *CharacterArray, n int) {
	a.n = n
	setU(a)
}

// IncludeArray merges two character arrays.
//
// It merges b into a, modifying a and returning a as the result.
// Merge means that the result has taxa that are in either a or b.
// If b.Len() > a.Len(), a will expand to the length of b.
func (a *CharacterArray) IncludeArray(b *CharacterArray) *CharacterArray {
	if n := b.Len(); n > a.n {
		setN(a, n)
	}
	a.bits.Or(a.bits, b.bits)
	return a
}

// IncludeTaxon sets the presence of the specified taxon in the reciever.
//
// The argument i corresponds to an index of the string representation.
// Values of i > a.Len() will expand a as needed.
//
// For convenience, the method returns the modified reciever.
func (a *CharacterArray) IncludeTaxon(i int) *CharacterArray {
	if i >= a.n {
		setN(a, i+1)
	}
	a.bits.SetBit(a.bits, i, 1)
	return a
}

// Len returns the total number of taxa (present or not) represented in the
// character array.
//
// See TaxaPresent for only the number of taxa present.
func (a *CharacterArray) Len() int {
	return a.n
}

// SetString sets the character array from a Go string.
//
// A '1' in the input string indicates the presence of the respective taxon.
// Any other byte indicates absence of the taxon.  If len(s) > a.Len(), a
// will be expanded as needed.  If len(s) < a.Len, taxa missing from the end
// of s are set to not present.
func (a *CharacterArray) SetString(s string) *CharacterArray {
	if n := len(s); n > a.n {
		setN(a, n)
	}
	a.bits.SetInt64(0)
	for bx := 0; bx < len(s); bx++ {
		if s[bx] == '1' {
			a.bits.SetBit(a.bits, bx, 1)
		}
	}
	return a
}

// String satisfies fmt.Stringer
func (a *CharacterArray) String() string {
	b := make([]byte, a.n)
	for i := range b {
		b[i] = '0' + byte(a.bits.Bit(i))
	}
	return string(b)
}

// TaxaPresent returns the number of taxa present in the character array.
//
// This is the number of 1's in the printed representation.
//
// See Len for the total number of taxa, 1's and 0's.
func (a *CharacterArray) TaxaPresent() (n int) {
	// yeah, there might be tricky ways to do this sub O(n).
	// it seems a premature optimization for now.
	for i := a.bits.BitLen(); i > 0; {
		i--
		if a.bits.Bit(i) == 1 {
			n++
		}
	}
	return
}

// Trivial returns true if the character array represents a trivial split.
// That is, A split with < 2 taxa present or with < 2 taxa missing.
func (a *CharacterArray) Trivial() bool {
	s := a.TaxaPresent()
	return s < 2 || a.n-s < 2
}

/* a kind of copy constructor. returns a new CharacterArray, equivalent to receiver
// but with bits inverted.
func (a *CharacterArray) Not() *CharacterArray {
	return &CharacterArray{a.n, a.universe, new(big.Int).AndNot(a.universe, a.bits)}
}*/

/*/ returns true if b is subset of a
func (a *CharacterArray) Subset(b *CharacterArray) bool {
	return new(big.Int).AndNot(b.bits, a.bits).BitLen() == 0
}*/

/* constructs new character array = a - b
func (a *CharacterArray) Subtract(b *CharacterArray) *CharacterArray {
	return &CharacterArray{a.n, a.universe, new(big.Int).AndNot(a.bits, b.bits)}
}*/

// CharacterTableFromStrings produces a phylogenic character table from
// a set of equal length symbol strings.
//
// At least four strings are required.  An error is returned for < 4 strings
// or for strings of unequal length.
func CharacterTableFromStrings(ss [][]byte) (r []*CharacterArray, err error) {
	if len(ss) < 4 {
		return nil, errors.New("not enough strings to characterize")
	}
	s0 := ss[0]
	sLen := len(ss[0])
	if sLen == 0 {
		return nil, errors.New("can't characterize empty strings")
	}
	var ca *CharacterArray
	for i := range s0 {
		if ca == nil {
			ca = NewCharacterArray(len(ss), nil, "1")
		} else {
			ca.SetString("1")
		}
		// oneSym is assigned immediately.  zeroSym == oneSym means that
		// a distinct symbol to correspond to character 0 has not been
		// encountered yet.
		oneSym := s0[i]   // symbol corresponding to character 1
		zeroSym := oneSym // symbol corresponding to character 0
	compareSym:
		for j := 1; j < len(ss); j++ {
			s := ss[j]
			if len(s) != sLen {
				return nil, errors.New("strings different lengths")
			}
			switch {
			case s[i] == oneSym:
				ca.IncludeTaxon(j)
			case zeroSym == oneSym:
				zeroSym = s[i]
			case s[i] == zeroSym:
			default:
				break compareSym
			}
		}
		if !ca.Trivial() {
			r = append(r, ca)
			ca = nil
		}
	}
	return
}

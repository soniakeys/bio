package bio

import (
	"testing"
)

func TestNewKMP(t *testing.T) {
	k := NewKMP([]byte("ababaa"))
	expected := []int{-1, 0, 0, 1, 2, 3, 1}
	for i, s := range k.b {
		if s != expected[i] {
			t.Fatal("got:\n", k, "\nexpected:\n", expected)
		}
	}
}

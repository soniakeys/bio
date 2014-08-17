package bio_test

import (
	"testing"

	"github.com/soniakeys/bio"
)

func TestNewKMP(t *testing.T) {
	k := bio.NewKMP([]byte("CAGTAAGCAGGGACTG"))
	expected := []int{0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 1, 0, 0}
	for i, s := range k {
		if s != expected[i] {
			t.Fatal("got:\n", k, "\nexpected:\n", expected)
		}
	}
}

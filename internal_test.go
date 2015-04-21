package bio

import "testing"

func TestLimbWeight(t *testing.T) {
	d := DistanceMatrix{
		{0, 13, 21, 22},
		{13, 0, 12, 13},
		{21, 12, 0, 13},
		{22, 13, 13, 0},
	}
	min, i, k := d.limbWeight(1)
	if min != 2 || i != 2 || k != 0 {
		t.Fatalf("got %d %d %d, want 2 2 0", min, i , k)
	}
}

func TestLimbWeightSubMatrix(t *testing.T) {
	d := DistanceMatrix{
		{0, 13, 21, 22},
		{13, 0, 12, 13},
		{21, 12, 0, 13},
		{22, 13, 13, 0},
	}
	min, i, k := d.limbWeightSubMatrix(3)
	if min != 7 || i != 1 || k != 2 {
		t.Fatalf("got %d %d %d, want 7 1 2", min, i , k)
	}
	min, i, k = d.limbWeightSubMatrix(2)
	if min != 10 || i != 0 || k != 1 {
		t.Fatalf("got %d %d %d, want 10 0 1", min, i , k)
	}
}

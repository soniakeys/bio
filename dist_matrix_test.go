package bio_test

import (
	"fmt"
	"testing"

	"github.com/soniakeys/bio"
)

func ExampleDistanceMatrix_LimbWeight() {
	d := bio.DistanceMatrix{
		{0, 13, 21, 22},
		{13, 0, 12, 13},
		{21, 12, 0, 13},
		{22, 13, 13, 0},
	}
	fmt.Println(d.LimbWeight(1))
	// Output:
	// 2 2 0
}

func ExampleDistanceMatrix_LimbWeightSubMatrix() {
	d := bio.DistanceMatrix{
		{0, 13, 21, 22},
		{13, 0, 12, 13},
		{21, 12, 0, 13},
		{22, 13, 13, 0},
	}
	fmt.Println(d.LimbWeightSubMatrix(3))
	fmt.Println(d.LimbWeightSubMatrix(2))
	// Output:
	// 7 1 2
	// 10 0 1
}

func ExampleDistanceMatrix_Valid() {
	d := bio.DistanceMatrix{
		{0, 13, 21, 22},
		{13, 0, 12, 13},
		{21, 12, 0, 13},
		{22, 13, 13, 0},
	}
	fmt.Println(d.Validate())
	// Output:
	// <nil>
}

func ExampleDistanceMatrix_Additive() {
	a := bio.DistanceMatrix{
		{0, 13, 21, 22},
		{13, 0, 12, 13},
		{21, 12, 0, 13},
		{22, 13, 13, 0},
	}
	na := bio.DistanceMatrix{
		{0, 3, 4, 3},
		{3, 0, 4, 5},
		{4, 4, 0, 2},
		{3, 5, 2, 0},
	}
	fmt.Println(a.Additive())
	fmt.Println(na.Additive())
	// Output:
	// true 0 0 0 0
	// false 3 1 0 2
}

func TestRandomAdditiveMatrix(t *testing.T) {
	for _, d := range []bio.DistanceMatrix{
		bio.RandomAdditiveMatrix(10),
		bio.RandomAdditiveMatrix(20),
		bio.RandomAdditiveMatrix(40),
	} {
		if ok, _, _, _, _ := d.Additive(); !ok {
			t.Fatal(len(d))
		}
	}
}

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

func ExampleDistanceMatrix_UPGMA() {
	d := bio.DistanceMatrix{
		{0, 20, 17, 11},
		{20, 0, 20, 13},
		{17, 20, 0, 10},
		{11, 13, 10, 0},
	}
	pl := d.UPGMA()
	fmt.Println("node  parent  weight     age  leaves")
	for i, u := range pl {
		fmt.Printf(">%3d     %3d  %6.3f  %6.3f     %3d\n",
			i, u.Parent, u.Weight, u.Age, u.NLeaves)
	}
	// Output:
	// node  parent  weight     age  leaves
	// >  0       5   7.000   0.000       1
	// >  1       6   8.833   0.000       1
	// >  2       4   5.000   0.000       1
	// >  3       4   5.000   0.000       1
	// >  4       5   2.000   5.000       2
	// >  5       6   1.833   7.000       3
	// >  6      -1     NaN   8.833       4
}

func ExampleDistanceMatrix_NeighborJoin() {
	d := bio.DistanceMatrix{
		{0, 23, 27, 20},
		{23, 0, 30, 28},
		{27, 30, 0, 30},
		{20, 28, 30, 0},
	}
	tree, wt := d.NeighborJoin()
	fmt.Println("n1  n2  weight")
	for n, to := range tree {
		for _, h := range to {
			fmt.Printf("%d  %2d   %6.3f\n", n, h.To, wt[h.Label])
		}
	}
	// Output:
	// n1  n2  weight
	// 0   5    8.000
	// 1   4   13.500
	// 2   4   16.500
	// 3   5   12.000
	// 4   5    2.000
	// 4   1   13.500
	// 4   2   16.500
	// 5   3   12.000
	// 5   0    8.000
	// 5   4    2.000
}

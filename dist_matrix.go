package bio

import (
	"errors"
	"fmt"
	"math/big"
	"math/rand"
)

// DistanceMatrix holds a distance matrix for phylogenic analyses.
//
// See DistanceMatrix.Valid for typical restrictions.
type DistanceMatrix [][]float64

// Valid validates a DistanceMarix.
//
// Conditions are:
//
//   * square:  len(d[i]) == len(d)
//   * symmetric:  d[i][j] == d[j][i]
//   * non-negative:  d[i][j] >= 0
//   * zero diagonal:  d[i][i] == 0
//   * triangle inequality:  d[i][j] + d[j][k] <= d[i][k]
//
// Valid returns nil if all conditions are met, otherwise an error citing
// a condition not met.
func (d DistanceMatrix) Validate() error {
	for i, di := range d {
		if len(di) != len(d) {
			return errors.New("not square")
		}
		for j, dij := range di {
			if dij < 0 {
				return fmt.Errorf("negative element: %g", dij)
			}
			if !(dij == d[j][i]) { // reversed test catches NaNs too.
				return errors.New("not symmetric")
			}
		}
		if !(di[i] == 0) {
			return errors.New("non-zero diagonal")
		}
	}
	for i, di := range d {
		for k, dk := range d[:i] {
			dik := di[k]
			for j, dij := range di {
				if dij+dk[j] < dik {
					return fmt.Errorf("triangle inequality not satisfied: "+
						"d[%d][%d] + d[%d][%d] < d[%d][%d]", i, j, j, k, i, k)
				}
			}
		}
	}
	return nil
}

// Additive tests if DistanceMatrix d is additive.
//
// Additive tests the four-point condition for all combinations of points.
// If a test fails, it returns ok = false and the four failing points.
func (d DistanceMatrix) Additive() (ok bool, i, j, k, l int) {
	for i, di := range d {
		for j, dj := range d[:i] {
			dij := di[j]
			for k, dk := range d[:j] {
				dik := di[k]
				djk := dj[k]
				for l, dil := range di {
					s1 := dij + dk[l]
					s2 := dik + dj[l]
					s3 := dil + djk
					switch { // swap so s1 is the value !=
					case s1 == s2:
						s1, s3 = s3, s1
					case s1 == s3:
						s1, s2 = s2, s1
					}
					if s1 > s2 {
						return false, i, j, k, l
					}
				}
			}
		}
	}
	ok = true
	return
}

// LimbWeight finds the weight of the edge to a leaf (a limb) in a phylogenic
// tree corresponding to DistanceMatrix d.
//
// Argument j is an index of d representing a leaf node.
//
// Returned is the weight wt of an edge that would connect j to a phylogenic
// tree, and two other leaf indexes, i and k, such that the edge connecting j
// to the tree would connect directly to the path between i and k.
func (d DistanceMatrix) LimbWeight(j int) (wt float64, i, k int) {
	// algorithm: pick k for convenience in iteration, iterate i over all
	// other indexes, accumulating the value of i that minimizes wt.
	var dj0, dj1, dk0, dk1 []float64
	if j == 0 {
		k = 1
		dj1 = d[j][2:]
		dk1 = d[k][2:]
	} else {
		k = j - 1
		dj0 = d[j][:k]
		dk0 = d[k][:k]
		dj1 = d[j][j+1:]
		dk1 = d[k][j+1:]
	}
	djk := d[j][k]
	var iMin2 int      // index i minimizing wtMin2
	var wtMin2 float64 // accumulated sum, twice min, return wtMin2/2
	if len(dj0) > 0 {
		iMin2 = len(dj0) - 1
		dj0 = dj0[:iMin2]
		dk0 = dk0[:iMin2]
	} else {
		iMin2 = 2
		dj1 = d[j][j+2:]
		dk1 = d[k][j+2:]
	}
	wtMin2 = d[j][iMin2] + djk - d[k][iMin2]
	for i0, dij := range dj0 {
		if wt2 := dij + djk - dk0[i]; wt2 < wtMin2 {
			wtMin2 = wt2
			iMin2 = i0
		}
	}
	for i1, dij := range dj1 {
		if wt2 := dij + djk - dk1[i1]; wt2 < wtMin2 {
			wtMin2 = wt2
			iMin2 = len(dj0) + 3 + i1
		}
	}
	return wtMin2 / 2, iMin2, k
}

// LimbWeightSubMatrix finds the weight of the edge to a leaf (a limb) in
// a phylogenic tree corresponding to a submatrix of d.
//
// Argument j is an index of d representing a leaf node.  The submatrix
// considered is d[:j+1][:j+1].  This variation of LimbWeight allows
// simpler and more efficient code.
//
// Return values are the same as for LimbWeight.
func (d DistanceMatrix) LimbWeightSubMatrix(j int) (wt float64, i, k int) {
	k = j - 1
	dj := d[j]
	dk := d[k]
	djk := dj[k]
	iMin2 := k - 1
	wtMin2 := dj[iMin2] + djk - dk[iMin2]
	for i, dij := range dj[:iMin2] {
		if wt := dij + djk - dk[i]; wt < wtMin2 {
			wtMin2 = wt
			iMin2 = i
		}
	}
	return wtMin2 / 2, iMin2, k
}

// A Half is an element in an adjacency list, a graph representation.
type Half struct {
	To     int
	Weight float64
}

// AdditiveTree constructs a phylogenic tree from an additive distance matrix.
//
// DistanceMatrix d must be additive.  Use provably additive matrices or
// use DistanceMatrix.Additive() to verify the additive property.
//
// Result is an unrooted tree, not necessarily binary, as an undirected graph.
// The structure is an adjacency list of Half structs.  The first len(d)
// nodes are the leaves represented by the distance matrix.  Internal nodes
// follow.
//
// Time complexity is O(n^2) in the number of leaves.
func (d DistanceMatrix) AdditiveTree() [][]Half {
	// interpretation of the presented recursive algorithm.
	// good things to try:  1: construct result as a parent list rather than
	// a child tree.  2: drop the recursion.  3. make tree always binary.
	t := make([][]Half, len(d)) // allocate leaves
	var ap func(int)
	ap = func(n int) {
		if n == 1 {
			wt := d[0][1]
			t[0] = []Half{{1, wt}}
			t[1] = []Half{{0, wt}}
			return
		}
		nLen, i, k := d.LimbWeightSubMatrix(n)
		x := d[i][n] - nLen
		ap(n - 1)
		// f() finds and returns connection node v.
		// method: df search to find i from k, find connection point on the
		// way out.
		// create connection node v if needed, return v if found, -1 if not.
		var vis big.Int
		var f func(n int) int
		f = func(n int) int {
			if n == i {
				return i
			}
			vis.SetBit(&vis, n, 1)
			for tx, to := range t[n] {
				if vis.Bit(to.To) == 1 {
					continue
				}
				p := f(to.To)
				switch {
				case p < 0: // not found yet
					continue
				case x == 0: // p is connection node
					return p
				case x < to.Weight: // new node at distance x from to.to
					v := len(t)
					y := to.Weight - x
					t[n][tx] = Half{v, y}
					for fx, from := range t[to.To] {
						if from.To == n {
							t[to.To][fx] = Half{v, x}
							break
						}
					}
					t = append(t, []Half{{n, y}, {to.To, x}})
					x = 0
					return v
				default: // continue back out
					x -= to.Weight
					return n
				}
			}
			return -1
		}
		vis = big.Int{}
		v := f(k)
		t[n] = []Half{{v, nLen}}
		t[v] = append(t[v], Half{n, nLen})
	}
	ap(len(d) - 1)
	return t
}

// A Half is an element in a parent list, an efficient representation for trees.
type FromHalf struct {
	From   int
	Weight float64
}

// Random binary tree, returned as a parent list.
//
// The binary tree for n leaves has n-2 internal nodes but the root node
// is not stored so len(parentList) == nLeaves+nLeaves-3.
// The leaves correspond to nodes 0:nLeaves and the unrepresented root
// corresponds to len(parentList).
//
// Example code:
//   n := 5
//   pl := bio.RandomBinaryTree(n)
//   fmt.Println(n+n-2, "nodes,")
//   fmt.Println(len(pl), "in parent list:")
//   for i, fh := range pl {
//       fmt.Printf("%d: parent %d weight %g\n", i, fh.From, fh.Weight)
//   }
//
// Output:
//
//   8 nodes,
//   7 in parent list:
//   0: parent 7 weight 62
//   1: parent 5 weight 52
//   2: parent 6 weight 24
//   3: parent 5 weight 84
//   4: parent 6 weight 98
//   5: parent 7 weight 82
//   6: parent 7 weight 84
//
// Diagram:
//
//          7
//        / | \
//       /  /  6
//      /  /   |\
//     /  5    | \
//    /  / \   |  \
//   0  1   3  2   4
//
func RandomBinaryTree(nLeaves int) (parentList []FromHalf) {
	// allocate space for whole tree except root
	parentList = make([]FromHalf, nLeaves+nLeaves-3)
	randWt := func() float64 {
		return float64(10 + rand.Intn(90))
	}
	// initial tree has three leaves and the internal root
	root := len(parentList)
	parentList[0] = FromHalf{root, randWt()}
	parentList[1] = FromHalf{root, randWt()}
	parentList[2] = FromHalf{root, randWt()}
	// now add remaining nodes.  on each iteration add one each internal
	// node and parent node.  pick random edge from parent list,
	// embed a new internal node in this edge, also connect a new leaf.
	// new edges are from new leaf to new internal node and from
	// new internal node to parent.
	for newLeaf := 3; newLeaf < nLeaves; newLeaf++ {
		i := nLeaves + newLeaf - 3     // new internal node
		l1 := rand.Intn(newLeaf*2 - 3) // (range is number of existing edges)
		if l1 >= newLeaf {
			l1 += nLeaves - newLeaf // skip to range of internal nodes
		}
		l2 := parentList[l1].From
		parentList[l1].From = i
		parentList[i] = FromHalf{l2, randWt()}      // new edge
		parentList[newLeaf] = FromHalf{i, randWt()} // new edge
	}
	return
}

// RAMatrix constructs a random additive distance matrix.
//
// Argument n is the size of the DistanceMatrix to reutrn.
func RandomAdditiveMatrix(n int) DistanceMatrix {
	pl := RandomBinaryTree(n)
	da := make([]struct { // distance annotation of parent list
		leng int     // path length
		dist float64 // distance to root
	}, len(pl))
	var f func(int) (int, float64)
	f = func(n int) (int, float64) {
		switch {
		case n == len(pl):
			return 1, 0
		case da[n].leng > 0:
			return da[n].leng, da[n].dist
		}
		leng, dist := f(pl[n].From)
		leng++
		dist += pl[n].Weight
		da[n].leng = leng
		da[n].dist = dist
		return leng, dist
	}
	for leaf := 0; leaf < n; leaf++ {
		f(leaf)
	}
	// distance between leaves in annotated parent list.
	lldist := func(l1, l2 int) (d float64) {
		// make l1 the leaf with the longer path
		if da[l1].leng < da[l2].leng {
			l1, l2 = l2, l1
		}
		// accumulate l1 distance to same tree height as l2
		len2 := da[l2].leng
		for da[l1].leng > len2 {
			d += pl[l1].Weight
			l1 = pl[l1].From
		}
		// accumulate d1, d2 until l1, l2 are the same node
		for l1 != l2 {
			d += pl[l1].Weight + pl[l2].Weight
			l1 = pl[l1].From
			l2 = pl[l2].From
		}
		return
	}
	// build distance matrix
	d := make(DistanceMatrix, n)
	for i := range d {
		di := make([]float64, n)
		for j := 0; j < i; j++ {
			di[j] = lldist(i, j)
			d[j][i] = di[j]
		}
		d[i] = di
	}
	return d
}

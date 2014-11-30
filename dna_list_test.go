package bio_test

import (
	"fmt"
	"math"

	"github.com/soniakeys/bio"
)

func ExampleKmers_Entropy() {
	k := bio.Kmers{
		bio.DNA8("GtCG"),
		bio.DNA8("GCtG"),
		bio.DNA8("GCCt"),
		bio.DNA8("cCCG"),
		bio.DNA8("GgCG"),
	}
	fmt.Printf("%.2f\n", k.Entropy())
	// Output:
	// 3.54
}

func ExampleKmers_Entropies() {
	k := bio.Kmers{
		bio.DNA8("GtCG"),
		bio.DNA8("GCtG"),
		bio.DNA8("GCCt"),
		bio.DNA8("cCCG"),
		bio.DNA8("GgCG"),
	}
	for i, e := range k.EntropyContributions() {
		fmt.Printf("Position %d: %.2f\n", i+1, e)
	}
	// Output:
	// Position 1: [0.00 0.46 0.00 0.26]
	// Position 2: [0.00 0.44 0.46 0.46]
	// Position 3: [0.00 0.26 0.46 0.00]
	// Position 4: [0.00 0.00 0.46 0.26]
}

func ExampleFracProfile_Entropy() {
	k := bio.Kmers{
		bio.DNA8("GtCG"),
		bio.DNA8("GCtG"),
		bio.DNA8("GCCt"),
		bio.DNA8("cCCG"),
		bio.DNA8("GgCG"),
	}
	fmt.Printf("%.2f\n", k.FracProfile().Entropy())
	// Output:
	// 3.54
}

func ExampleFracProfile_CrossEntropy() {
	d := bio.DNA8List{
		bio.DNA8("taaaaGtCGa"),
		bio.DNA8("acGCtGaaaa"),
		bio.DNA8("aaaaGCCtat"),
		bio.DNA8("acCCGaataa"),
		bio.DNA8("agaaaaGgCG"),
	}
	b := d.BaseFreq()
	for i, f := range b {
		b[i] = math.Log2(f)
	}
	k := bio.Kmers{
		bio.DNA8("GtCG"),
		bio.DNA8("GCtG"),
		bio.DNA8("GCCt"),
		bio.DNA8("cCCG"),
		bio.DNA8("GgCG"),
	}
	fmt.Printf("%.2f\n", k.FracProfile().CrossEntropy(b))
	// Output:
	// 9.97
}

func ExampleFracProfile_RelativeEntropy() {
	d := bio.DNA8List{
		bio.DNA8("taaaaGtCGa"),
		bio.DNA8("acGCtGaaaa"),
		bio.DNA8("aaaaGCCtat"),
		bio.DNA8("acCCGaataa"),
		bio.DNA8("agaaaaGgCG"),
	}
	b := d.BaseFreq()
	k := bio.Kmers{
		bio.DNA8("GtCG"),
		bio.DNA8("GCtG"),
		bio.DNA8("GCCt"),
		bio.DNA8("cCCG"),
		bio.DNA8("GgCG"),
	}
	fmt.Printf("%.2f\n", k.FracProfile().RelativeEntropy(b))
	// Output:
	// 6.44
}

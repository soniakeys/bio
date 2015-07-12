package bio_test

import (
	"fmt"
	"math"
	"sort"

	"github.com/soniakeys/bio"
)

func ExampleDNA8List() {
	l := bio.DNA8List{
		bio.DNA8("aAcCtTgG"),
		bio.DNA8("a"),
		bio.DNA8("actg"),
		bio.DNA8("aCtg"),
		bio.DNA8("acTg"),
		bio.DNA8("acgt"),
	}
	sort.Sort(l)
	for _, seq := range l {
		fmt.Println(seq)
	}
	// Output:
	// a
	// actg
	// acTg
	// acgt
	// aCtg
	// aAcCtTgG
}

func ExampleKmers() {
	k := bio.Kmers{
		bio.DNA8("actg"),
		bio.DNA8("aCtg"),
		bio.DNA8("acTg"),
		bio.DNA8("acgt"),
	}
	l := bio.DNA8List(k)
	sort.Sort(l)
	for _, seq := range l {
		fmt.Println(seq)
	}
	// Output:
	// actg
	// acTg
	// acgt
	// aCtg
}

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

func ExampleMedianMotifs() {
	d := bio.DNA8List{
		bio.DNA8("AAATTGACGCAT"),
		bio.DNA8("GACGACCACGTT"),
		bio.DNA8("CGTCAGCGCCTG"),
		bio.DNA8("GCTGAGCACCGG"),
		bio.DNA8("AGTTCGGGACAG"),
	}
	fmt.Println(d.MedianMotifs(3))
	// Output:
	// [GAC] 2
}

func ExampleDNA8List_GreedyMotifSearch() {
	d := bio.DNA8List{
		bio.DNA8("GGCGTTCAGGCA"),
		bio.DNA8("AAGAATCAGTCA"),
		bio.DNA8("CAAGGAGTTCGC"),
		bio.DNA8("CACGTCAATCAC"),
		bio.DNA8("CAATAATATTCG"),
	}
	fmt.Println(d.GreedyMotifSearch(3))
	// Output:
	// [[TTC ATC TTC ATC TTC] [TCA TCA TCG TCA TCG] [CAG CAG CAA CAA CAA]] 2
}

package bio_test

import (
	"fmt"
	"math"
	"sort"
	"testing"

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

func ExampleDNA8List_PlantedMotifsPMS1() {
	d := bio.DNA8List{
		bio.DNA8("AAATTGACGCAT"),
		bio.DNA8("GACGACCACGTT"),
		bio.DNA8("CGTCAGCGCCTG"),
		bio.DNA8("GCTGAGCACCGG"),
		bio.DNA8("AGTTCGGGACAG"),
	}
	for m := range d.PlantedMotifsPMS1(4, 1) {
		fmt.Println(m)
	}
	// Output:
	// GCAC
}

func ExampleDNA8List_PlantedMotifs() {
	d := bio.DNA8List{
		bio.DNA8("AAATTGACGCAT"),
		bio.DNA8("GACGACCACGTT"),
		bio.DNA8("CGTCAGCGCCTG"),
		bio.DNA8("GCTGAGCACCGG"),
		bio.DNA8("AGTTCGGGACAG"),
	}
	for _, m := range d.PlantedMotifs(4, 1) {
		fmt.Println(m)
	}
	// Output:
	// GCAC
}

func ExampleDNA8List_MedianMotifs() {
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

func ExampleDNA8List_RandomMotifSearch() {
	d := bio.DNA8List{
		bio.DNA8("CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA"),
		bio.DNA8("GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG"),
		bio.DNA8("TAGTACCGAGACCGAAAGAAGTATACAGGCGT"),
		bio.DNA8("TAGATCAAGTTTCAGGTGCACGTCGGTGAACC"),
		bio.DNA8("AATCCACCAGCTCCACGTGCAATGTTGGCCTA"),
	}
	fmt.Println(d.RandomMotifSearch(8, 1000))
	// Output:
	// [[TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG]] 9
}

func ExampleDNA8List_GibbsMotifSearch() {
	d := bio.DNA8List{
		bio.DNA8("CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA"),
		bio.DNA8("GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG"),
		bio.DNA8("TAGTACCGAGACCGAAAGAAGTATACAGGCGT"),
		bio.DNA8("TAGATCAAGTTTCAGGTGCACGTCGGTGAACC"),
		bio.DNA8("AATCCACCAGCTCCACGTGCAATGTTGGCCTA"),
	}
	ks, h := d.GibbsMotifSearch(8, 100, 20)
	fmt.Println(h)
	// (sort text for predictable output)
	s := make([]string, len(ks))
	for i, k := range ks {
		s[i] = fmt.Sprint(k)
	}
	sort.Strings(s)
	for _, k := range s {
		fmt.Println(k)
	}
	// Output:
	// 9
	// [AACGGCCA AAGTGCCA TAGTACCG AAGTTTCA ACGTGCAA]
	// [TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG]
}

var m10 = bio.DNA8List{
	bio.DNA8("TGATGATAACGTGACGGGACTCAGCGGCGATGAAGGATGAGT"),
	bio.DNA8("CAGCGACAGACAATTTCAATAATATCCGCGGTAAGCGGCGTA"),
	bio.DNA8("TGCAGAGGTTGGTAACGCCGGCGACTCGGAGAGCTTTTCGCT"),
	bio.DNA8("TTTGTCATGAACTCAGATACCATAGAGCACCGGCGAGACTCA"),
	bio.DNA8("ACTGGGACTTCACATTAGGTTGAACCGCGAGCCAGGTGGGTG"),
	bio.DNA8("TTGCGGACGGGATACTCAATAACTAAGGTAGTTCAGCTGCGA"),
	bio.DNA8("TGGGAGGACACACATTTTCTTACCTCTTCCCAGCGAGATGGC"),
	bio.DNA8("GAAAAAACCTATAAAGTCCACTCTTTGCGGCGGCGAGCCATA"),
	bio.DNA8("CCACGTCCGTTACTCCGTCGCCGTCAGCGATAATGGGATGAG"),
	bio.DNA8("CCAAAGCTGCGAAATAACCATACTCTGCTCAGGAGCCCGATG"),
}

func BenchmarkMedian(b *testing.B) {
	for i := 0; i < b.N; i++ {
		m10.MedianMotifs(6)
	}
}

func BenchmarkMedianB(b *testing.B) {
	for i := 0; i < b.N; i++ {
		m10.MedianMotifsB(6)
	}
}

/* 4 1
var plants = bio.DNA8List{
	bio.DNA8("AAATTGACGCAT"),
	bio.DNA8("GACGACCACGTT"),
	bio.DNA8("CGTCAGCGCCTG"),
	bio.DNA8("GCTGAGCACCGG"),
	bio.DNA8("AGTTCGGGACAG"),
}
*/

var plants = bio.DNA8List{
	bio.DNA8("TCTGAGCTTGCGTTATTTTTAGACC"),
	bio.DNA8("GTTTGACGGGAACCCGACGCCTATA"),
	bio.DNA8("TTTTAGATTTCCTCAGTCCACTATA"),
	bio.DNA8("CTTACAATTTCGTTATTTATCTAAT"),
	bio.DNA8("CAGTAGGAATAGCCACTTTGTTGTA"),
	bio.DNA8("AAATCCATTAAGGAAAGACGACCGT"),
}

func BenchmarkPlantedMotifs(b *testing.B) {
	for i := 0; i < b.N; i++ {
		plants.PlantedMotifs(5, 2)
	}
}

func BenchmarkPlantedMotifs1(b *testing.B) {
	for i := 0; i < b.N; i++ {
		plants.PlantedMotifs1(5, 2)
	}
}

func BenchmarkPlantedMotifs2(b *testing.B) {
	for i := 0; i < b.N; i++ {
		plants.PlantedMotifs2(5, 2)
	}
}

func BenchmarkPlantedMotifsPMS1(b *testing.B) {
	for i := 0; i < b.N; i++ {
		plants.PlantedMotifsPMS1(5, 2)
	}
}

func ExampleDNA8List_KCompositionDistMat() {
	l := bio.DNA8List{
		bio.DNA8("ATATATAG"),
		bio.DNA8("ATATATA"),
		bio.DNA8("GATATA"),
	}
	d := l.KCompositionDistMat(3)
	for i, s := range l {
		fmt.Printf("%-8s %.0f\n", s, d[i])
	}
	// Output:
	// ATATATAG [0 1 4]
	// ATATATA  [1 0 3]
	// GATATA   [4 3 0]
}

func ExampleDNA8List_DistanceMatrix() {
	l := bio.DNA8List{
		bio.DNA8("TTTCCATTTA"),
		bio.DNA8("GATTCATTTC"),
		bio.DNA8("TTTCCATTTT"),
		bio.DNA8("GTTCCATTTA"),
	}
	f := func(a, b bio.DNA8) float64 { return float64(a.Hamming(b)) }
	m := l.DistanceMatrix(f)
	for _, r := range m {
		fmt.Printf("%2.0f\n", r)
	}
	// Output:
	// [ 0  4  1  1]
	// [ 4  0  4  3]
	// [ 1  4  0  2]
	// [ 1  3  2  0]
}

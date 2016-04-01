package bio_test

import (
	"fmt"
	"log"
	"sort"

	"github.com/soniakeys/bio"
)

func ExampleReadPairList_ReadBreak() {
	ps := bio.ReadPairList{2, []bio.KmerPair{
		{"ABCD", "GHIJ"},
		{"CDEF", "IJKL"},
	}}
	f, err := ps.ReadBreak(2)
	if err != nil {
		log.Fatal(err)
	}
	var s []string
	for p, n := range f.Freq {
		s = append(s, fmt.Sprint(p, " ", n))
	}
	sort.Strings(s)
	fmt.Println("d:", f.D)
	for _, l := range s {
		fmt.Println(l)
	}
	// Output:
	// d: 4
	// {AB GH} 1
	// {BC HI} 1
	// {CD IJ} 2
	// {DE JK} 1
	// {EF KL} 1
}

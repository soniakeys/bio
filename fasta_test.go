package bio_test

import (
	"bytes"
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleReadFasta() {
	r := bytes.NewBufferString(`>sequence header
AGACCA
TACCA`)
	s, err := bio.ReadFasta(r)
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Println(s[0].Header)
	fmt.Println(bio.DNA(s[0].Seq))
	// Output:
	// >sequence header
	// AGACCATACCA
}

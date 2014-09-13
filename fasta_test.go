package bio_test

import (
	"bytes"
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleReadFASTA() {
	r := bytes.NewBufferString(`>sequenceID more stuff
AGACCA
TACCA`)
	s, err := bio.ReadFASTA(r)
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Println(len(s), "sequence")
	f := s[0]
	fmt.Println(f.ID())
	fmt.Println(f.Header)
	fmt.Println(bio.DNA(s[0].Seq))
	// Output:
	// 1 sequence
	// sequenceID
	// >sequenceID more stuff
	// AGACCATACCA
}

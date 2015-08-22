package bio_test

import (
	"bytes"
	"fmt"

	"github.com/soniakeys/bio"
)

func ExampleReadFASTA() {
	r := bytes.NewBufferString(`>db|123|abc example sequence
AGACCA
TACCA`)
	s, err := bio.ReadFASTA(r)
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Println(len(s), "sequence:")
	f := s[0]
	fmt.Println("Header:", f.Header)
	fmt.Println("Seq:   ", f.Seq)
	// Output:
	// 1 sequence:
	// Header: >db|123|abc example sequence
	// Seq:    AGACCATACCA
}

func ExampleFASTASeq() {
	r := bytes.NewBufferString(`>db|123|abc example sequence
AGACCA
TACCA`)
	s, err := bio.ReadFASTA(r)
	if err != nil {
		fmt.Println(err)
		return
	}
	f := s[0]
	fmt.Println("ID:  ", f.ID())
	fmt.Println("Desc:", f.Desc())
	// Output:
	// ID:   db|123|abc
	// Desc: example sequence
}

func ExampleFASTAReader_ReadSeq() {
	b := bytes.NewBufferString(`>db|123|abc example sequence
AGACCA
TACCA`)
	r := bio.NewFASTAReader(b)
	s, err := r.ReadSeq()
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Println("Header:", s.Header)
	fmt.Println("Seq:   ", s.Seq)
	_, err = r.ReadSeq()
	fmt.Println(err)
	// Output:
	// Header: >db|123|abc example sequence
	// Seq:    AGACCATACCA
	// EOF
}

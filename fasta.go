package bio

import (
	"bytes"
	"errors"
	"io"
	"io/ioutil"
)

// FastaSeq is a sequence annotated with a header.
type FastaSeq struct {
	Header string
	Seq    []byte
}

// ReadFasta reads FASTA format from an io.Reader
//
// Super minimal whole-file reader.  It allows blank lines,
// it trims leading and trailing white space, and it
// allows both lf and crlf line endings.  No other frills.
// Headers are returned unparsed.  The leading > isn't even removed.
//
// The first non-blank line must be a header.  A FastaSeq is returned
// for every header line, even if there is no data following the header.
func ReadFasta(r io.Reader) (seq []FastaSeq, err error) {
	d, err := ioutil.ReadAll(r)
	if err != nil {
		return nil, err
	}
	var last *FastaSeq
	foundHeader := false
	for _, l := range bytes.Split(d, []byte{'\n'}) {
		l = bytes.TrimSpace(l)
		switch {
		case len(l) == 0:
			continue // skip blank line
		case l[0] == '>':
			foundHeader = true // new seq
			seq = append(seq, FastaSeq{string(l), nil})
			last = &seq[len(seq)-1]
			continue
		case !foundHeader:
			return nil, errors.New("No header")
		}
		last.Seq = append(last.Seq, l...) // store seq data
	}
	return seq, nil
}

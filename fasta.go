package bio

import (
	"bufio"
	"bytes"
	"errors"
	"io"
	"io/ioutil"
	"os"
	"strings"
)

// This file currently has two separate readers.  They both return FASTASeq
// results at this point, But need further integration.  (TODO)

// FASTASeq is a sequence annotated with a header.
type FASTASeq struct {
	Header string // complete header line, including '>'
	Seq    []byte
}

// ID extracts the sequence identifier from the header.
func (f FASTASeq) ID() string {
	if f.Header == "" {
		return ""
	}
	id := f.Header[1:]
	if sp := strings.IndexByte(id, ' '); sp > 0 {
		return id[:sp]
	}
	return id
}

// ReadFASTA reads FASTA format from an io.Reader
//
// Super minimal whole-file reader.  It allows blank lines,
// it trims leading and trailing white space, and it
// allows both lf and crlf line endings.  No other frills.
// Headers are stored unparsed.
//
// The first non-blank line must be a header.  A FASTASeq is returned
// for every header line, even if there is no data following the header.
func ReadFASTA(r io.Reader) (seq []FASTASeq, err error) {
	d, err := ioutil.ReadAll(r)
	if err != nil {
		return nil, err
	}
	var last *FASTASeq
	foundHeader := false
	for _, l := range bytes.Split(d, []byte{'\n'}) {
		l = bytes.TrimSpace(l)
		switch {
		case len(l) == 0:
			continue // skip blank line
		case l[0] == '>':
			foundHeader = true // new seq
			seq = append(seq, FASTASeq{string(l), nil})
			last = &seq[len(seq)-1]
			continue
		case !foundHeader:
			return nil, errors.New("No header")
		}
		last.Seq = append(last.Seq, l...) // store seq data
	}
	return seq, nil
}

// second reader follows

// FASTAReader type for representing a FASTA stream.
type FASTAReader struct {
	r          *bufio.Reader
	nextHeader string
}

// NewFASTAReader constructs and returns a FASTAReader around a bufio.Reader.
func NewFASTAReader(r *bufio.Reader) FASTAReader {
	return FASTAReader{r: r}
}

// ReadSequence returns a single sequence on each call.
//
// A successful read is indicated by err = nil for all sequences, including
// the last.  Subsequent calls return err = io.EOF.
// Other error values indicate problems.
func (r *FASTAReader) ReadSequence() (FASTASeq, error) {
	f := FASTASeq{r.nextHeader, nil}
	r.nextHeader = ""
	var line []byte
	var isPre bool
	var err error
read:
	for {
		line, isPre, err = r.r.ReadLine()
		switch {
		case err != nil:
			if err == io.EOF && f.Header > "" {
				err = nil
			}
			break read
		case isPre:
			err = errors.New("FASTA line too long")
			break read
		case len(line) == 0: // ignore blank lines
			continue
		case line[0] == '>': // info line
			if f.Header > "" {
				r.nextHeader = string(line)
				break read
			}
			f.Header = string(line)
		case f.Header > "":
			f.Seq = append(f.Seq, line...)
		default:
			// ignore inital data without Header?
		}
	}
	return f, err
}

// ReadFASTAFile is a high level function that reads an entire FASTA file
// into memory.
func ReadFASTAFile(path string) ([]FASTASeq, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	var list []FASTASeq
	for r := NewFASTAReader(bufio.NewReader(f)); ; {
		s, err := r.ReadSequence()
		if err != nil {
			break
		}
		list = append(list, s)
	}
	if err != io.EOF {
		return nil, err
	}
	return list, nil
}

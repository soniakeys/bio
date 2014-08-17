package bio

import (
	"bufio"
	"bytes"
	"errors"
	"io"
	"io/ioutil"
	"os"
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

// FASTAReader type for representing a FASTA stream.
type FASTAReader struct {
	r      *bufio.Reader
	nextID string
}

// NewFASTAReader constructs and returns a FASTAReader around a bufio.Reader.
func NewFASTAReader(r *bufio.Reader) FASTAReader {
	return FASTAReader{r: r}
}

// ReadSequence returns a single sequence on each call.
//
// The entire header line following '>' is returned as seqID.
// A successful read is indicated by err = nil for all sequences, including
// the last.  Subsequent calls return err = io.EOF.
// Other error values indicate problems.
func (r *FASTAReader) ReadSequence() (seqID string, seq []byte, err error) {
	seqID = r.nextID
	r.nextID = ""
	var line []byte
	var isPre bool
read:
	for {
		line, isPre, err = r.r.ReadLine()
		switch {
		case err != nil:
			if err == io.EOF && seqID > "" {
				err = nil
			}
			break read
		case isPre:
			err = errors.New("FASTA line too long")
			break read
		case len(line) == 0: // ignore blank lines
			continue
		case line[0] == '>': // info line
			if seqID > "" {
				r.nextID = string(line[1:])
				break read
			}
			seqID = string(line[1:])
		case seqID > "":
			seq = append(seq, line...)
		default:
			// ignore string data without seqID
		}
	}
	return
}

// FASTASeq is the type returned by ReadFASTA.
type FASTASeq struct {
	SeqID string
	Seq   []byte
}

// ReadFASTA is a high level function that reads an entire FASTA file into
// memory.
func ReadFASTA(path string) (list []FASTASeq, err error) {
	var f *os.File
	if f, err = os.Open(path); err != nil {
		return
	}
	var s FASTASeq
	for r := NewFASTAReader(bufio.NewReader(f)); ; {
		s.SeqID, s.Seq, err = r.ReadSequence()
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

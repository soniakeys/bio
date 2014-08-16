package bio_test

import (
	"testing"

	"github.com/soniakeys/bio"
)

func TestBlosum62(t *testing.T) {
	if bio.Blosum62 == nil {
		t.Fatal("Blosum62 not initialized")
	}
	if score := bio.Blosum62.LinearGap(
		bio.AA("PLEASANTLY"),
		bio.AA("-MEA--N-LY"), 5); score != 8 {
		t.Fatal("Blosum62 score:", score)
	}
}

func TestPAM250(t *testing.T) {
	if bio.PAM250 == nil {
		t.Fatal("PAM250 not initialized")
	}
	if score := bio.PAM250.LinearGap(
		bio.AA("LYPRTEINSTRIN"),
		bio.AA("LY---EINSTEIN"), 5); score != 23 {
		t.Fatal("PAM250 score:", score)
	}
}

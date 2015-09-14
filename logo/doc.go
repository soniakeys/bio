// Logo generates a "motif logo," a graphical representation of a profile.
//
// Example program that serves a logo so it can be rendered with a browser:
/*
     package main

     import (
        "log"
        "net/http"

        "github.com/soniakeys/bio"
	    "github.com/soniakeys/bio/logo"
     )

     func main() {
        k := bio.Kmers{
           bio.DNA8("TCGGGGgTTTtt"),
           bio.DNA8("cCGGtGAcTTaC"),
           bio.DNA8("aCGGGGATTTtC"),
           bio.DNA8("TtGGGGAcTTtt"),
           bio.DNA8("aaGGGGAcTTCC"),
        }
        fp := k.FracProfile()

        // generate svg, setting options to tighten up margins a bit
        lg := logo.Motif(fp, logo.Margins(5), logo.BotMargin(20))

        // serve to localhost:2003/logo
        http.HandleFunc("/logo", func(w http.ResponseWriter, req *http.Request) {
           w.Header().Set("Content-Type", "image/svg+xml")
           w.Write(lg)
        })
        if err := http.ListenAndServe(":2003", nil); err != nil {
           log.Fatal("ListenAndServe:", err)
        }
     }
*/
// Most functions in the package just set options.  See doc for Motif()
// as a starting point.
package logo

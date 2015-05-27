package main

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func stripTags(path string, tags []string) {
	f, err := os.Open(path)
	check(err)
	b, err := bam.NewReader(f, 3)
	check(err)

	out := os.Stdout
	check(err)

	w, err := bam.NewWriter(out, b.Header(), 2)
	check(err)

	btags := make([][]byte, len(tags))
	for i := range tags {
		btags[i] = []byte(tags[i])
	}

	for {
		rec, err := b.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			check(err)
		}
		newAux := make(sam.AuxFields, len(rec.AuxFields))
		copy(newAux, rec.AuxFields)
		for _, btag := range btags {
			for i, aux := range newAux {
				if bytes.Compare(aux[:2], btag) == 0 {
					copy(newAux[i:], newAux[i+1:])
					newAux[len(newAux)-1] = nil
					newAux = newAux[:len(newAux)-1]
					break
				}
			}

		}
		rec.AuxFields = newAux
		e := w.Write(rec)
		check(e)

	}
	w.Close()
}

func main() {

	bam := flag.String("bam", "", "path to bam to strip tags")
	flag.Parse()
	tags := flag.Args()
	if len(tags) == 0 || *bam == "" {
		fmt.Printf("send in names of tags to strip and path to bam file\n")
		fmt.Printf("e.g. strip-tags -bam some.bam XS AS MC")
		flag.PrintDefaults()
		return
	}
	stripTags(*bam, tags)
}

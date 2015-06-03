// call somatic variants with multiple tumor samples associated with 1 normal. This uses
// the genotype likelihood approach from bcbio and speedseq but allows for multiple tumor
// samples. Anything that meets the criteria will have a "PASS" or "." FILTER and a list
// of tumor samples with the somatic variant in INFO["SOMATIC"]
// Copied from bcbio and speedseq.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strings"

	"github.com/brentp/vcfgo"
	"github.com/brentp/xopen"
)

func getTumorLOD(gls []float64) float64 {
	lod := -1.0
	if len(gls) < 3 {
		return lod
	}
	for _, gl := range gls[1:] {
		t := gl - gls[0]
		if t > lod {
			lod = t
		}
	}
	return lod
}

func getNormalLOD(gls []float64, thresh float64) float64 {
	lod := 1e7
	if len(gls) < 3 {
		return thresh
	}
	for _, gl := range gls[1:] {
		t := gls[0] - gl
		if t < lod {
			lod = t
		}
	}
	return lod
}

func somaticLOD(normalGLs []float64, tumorGLs []float64, thresh float64) bool {
	tumorLOD := getTumorLOD(tumorGLs)
	normalLOD := getNormalLOD(normalGLs, thresh)
	return tumorLOD >= thresh && normalLOD >= thresh
}

func getFreq(refCount int, altCounts []int) float64 {
	ac := 0
	for _, a := range altCounts {
		ac += a
	}
	return float64(ac) / float64(refCount+ac)
}

var THRESH_RATIO = 2.7

func somaticFreqs(normal *vcfgo.SampleGenotype, tumor *vcfgo.SampleGenotype) bool {
	var nFreq, tFreq float64
	nrd, err := normal.RefDepth()
	if err != nil {
		nFreq = 0.0
	} else {
		nads, err := normal.AltDepths()
		if err != nil {
			nFreq = 0.0
		} else {
			nFreq = getFreq(nrd, nads)
		}
	}
	trd, err := tumor.RefDepth()
	if err != nil {
		tFreq = 0.0
	} else {
		tads, err := tumor.AltDepths()
		if err != nil {
			tFreq = 0.0
		} else {
			tFreq = getFreq(trd, tads)
		}
	}
	return nFreq <= 0.001 || nFreq <= tFreq/THRESH_RATIO

}

func Somatics(v *vcfgo.Variant, normalIdx int) []string {
	thresh := 3.5

	normal := v.Samples[normalIdx]
	somatics := make([]string, 0)
	for i, tumor := range v.Samples {
		if i == normalIdx {
			continue
		}
		if !somaticFreqs(normal, tumor) {
			continue
		}
		if somaticLOD(normal.GL, tumor.GL, thresh) {
			somatics = append(somatics, v.Header.SampleNames[i])
		}
	}
	return somatics
}

func check(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

func main() {
	index := flag.Int("index", 0, "0-based index of the normal sample")
	flag.Parse()
	vcfs := flag.Args()
	if len(vcfs) != 1 {
		fmt.Printf("---------------- call somatic variants present in any tumor sample ----------------\n")
		fmt.Printf("----- uses the method from bcbio and speedseq but for multiple tumor samples ------\n")
		fmt.Printf("%s -index 0 normal-and-tumors.vcf.gz\n", os.Args[0])
		flag.PrintDefaults()
		os.Exit(1)
	}

	fhr, err := xopen.Ropen(vcfs[0])
	check(err)

	rdr, err := vcfgo.NewReader(fhr, false)
	check(err)

	fhw, err := xopen.Wopen("-")
	check(err)

	hdr := rdr.Header
	hdr.Filters["NOT_SOMATIC"] = "not somatic between normal and any tumor"

	hdr.Infos["SOMATIC"] = &vcfgo.Info{Id: "SOMATIC", Number: "1", Type: "String", Description: "Tumor samples with somatic event"}

	wtr, err := vcfgo.NewWriter(fhw, hdr)
	check(err)

	log.Printf("using %s as the normal sample\n", hdr.SampleNames[*index])

	for v := rdr.Read(); v != nil; v = rdr.Read() {

		somatics := Somatics(v, *index)
		if len(somatics) > 0 {
			v.Info.Add("SOMATIC", strings.Join(somatics, "|"))
		} else {
			if v.Filter == "." || v.Filter == "PASS" {
				v.Filter = ""
			}
			if v.Filter != "" {
				v.Filter += ";"
			}
			v.Filter += "NOT_SOMATIC"
		}
		wtr.WriteVariant(v)
	}
	fhw.Close()
}

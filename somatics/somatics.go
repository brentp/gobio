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
	normalLOD := getNormalLOD(normalGLs, thresh)
	tumorLOD := getTumorLOD(tumorGLs)
	return tumorLOD >= thresh && normalLOD >= thresh
}

func getFreq(refCount int, altCounts []int) float64 {
	ac := 0
	for _, a := range altCounts {
		ac += a
	}
	return float64(ac) / float64(refCount+ac)
}

func somaticFreqs(normal *vcfgo.SampleGenotype, tumor *vcfgo.SampleGenotype, freqRatio float64) bool {
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
	return nFreq <= 0.001 || nFreq <= tFreq/freqRatio

}

func Somatics(v *vcfgo.Variant, normalIdx int, thresh float64, freqRatio float64, skipMissing bool) []string {

	normal := v.Samples[normalIdx]
	somatics := make([]string, 0)
	// if we don't have GLs for normal, dont call a somatic.
	if skipMissing && len(normal.GL) == 0 {
		return somatics
	}
	for i, tumor := range v.Samples {
		if i == normalIdx {
			continue
		}
		if !somaticFreqs(normal, tumor, freqRatio) {
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

func sampleIndex(samples []string, sample string) int {
	for p, s := range samples {
		if s == sample {
			return p
		}
	}
	return -1
}

func main() {
	normalStr := flag.String("normal", "", "sample name of the normal sample. Default is the first sample")
	threshold := flag.Float64("threshold", 3.5, "threshold for difference in GLs to call a somatic variant")
	freqRatio := flag.Float64("freq-ratio", 2.7, "frequency in the tumor must be at least this many times the frequency in the normal")
	onlySomatic := flag.Bool("only-somatic", false, "print only the PASSing somatic variants (default is to set a flag and print all variants")
	skipMissing := flag.Bool("skip-missing-normals", false, "skip variants with missing normal. default is flag these as somatic.")
	flag.Parse()
	vcfs := flag.Args()
	if len(vcfs) != 1 {
		fmt.Printf("-----------------------------------------------------------------------------------\n")
		fmt.Printf("---------------- call somatic variants present in any tumor sample ----------------\n")
		fmt.Printf("----- uses the method from bcbio and speedseq but for multiple tumor samples ------\n")
		fmt.Printf("-----------------------------------------------------------------------------------\n")
		fmt.Printf("%s normal-and-tumors.vcf.gz\n", os.Args[0])
		flag.PrintDefaults()
		os.Exit(1)
	}

	/*
		f, err := os.Create("cpu.pprof")
		if err != nil {
			panic(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	*/

	fhr, err := xopen.Ropen(vcfs[0])
	check(err)

	rdr, err := vcfgo.NewReader(fhr, false)
	check(err)

	normal := 0
	if *normalStr != "" {
		normal = sampleIndex(rdr.Header.SampleNames, *normalStr)
		if normal == -1 {
			log.Fatalf("sample %s: not found in vcf (%v)\n", *normalStr, rdr.Header.SampleNames)
		}
	}

	fhw, err := xopen.Wopen("-")
	check(err)

	hdr := rdr.Header
	hdr.Filters["NOT_SOMATIC"] = "not somatic between normal and any tumor"

	hdr.Infos["SOMATIC"] = &vcfgo.Info{Id: "SOMATIC", Number: "1", Type: "String", Description: "Tumor samples with somatic event"}

	wtr, err := vcfgo.NewWriter(fhw, hdr)
	check(err)

	log.Printf("using %s as the normal sample\n", hdr.SampleNames[normal])

	j := 0
	for v := rdr.Read(); v != nil; v = rdr.Read() {
		j++
		if *onlySomatic && !(v.Filter == "." || v.Filter == "PASS") {
			continue
		}
		somatics := Somatics(v, normal, *threshold, *freqRatio, *skipMissing)
		if len(somatics) > 0 {
			v.Info.Add("SOMATIC", strings.Join(somatics, "|"))
		} else {
			if *onlySomatic {
				continue
			}
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
	log.Println("VCF warnings:", rdr.Error())
	log.Printf("evaluated %d variants\n", j)
	fhw.Close()
}

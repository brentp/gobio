// chunkbam accepts a bam file and the number of regions to split it into
// and outputs BED of intervals split evenly by number of reads. This is
// useful for parallelization. It can process a 7GB RNA-Seq file (note that this
// will have very uneven coverage) in about 1 minute, 25 seconds of user time with
// 20 CPUs. It starts outputting rows ready for parallelization in < 5 seconds.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"time"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

func chunkify(path string, nChunks int) chan string {
	f, err := os.Open(path)
	if err != nil {
		log.Fatal(err)
	}
	b, err := bam.NewReader(f, 0)
	if err != nil {
		log.Fatal("error opening bam:", err)
	}
	ifr, err := os.Open(path + ".bai")
	if err != nil {
		log.Println("bam must be indexed")
		log.Fatal(err)
	}

	bai, err := bam.ReadIndex(ifr)
	log.Println("read index")
	if err != nil {
		log.Fatal(err)
	}

	refs := b.Header().Refs()
	totalMapped := uint64(0)
	for _, ref := range refs {
		s, ok := bai.ReferenceStats(ref.ID())
		if !ok {
			continue
		}
		totalMapped += s.Mapped
	}
	readsPerChunk := int(float64(totalMapped) / float64(nChunks))

	chunks := make(map[string]chan string, 8)
	for _, ref := range refs {
		_, ok := bai.ReferenceStats(ref.ID())
		if !ok {
			continue
		}
		chunks[ref.Name()] = make(chan string, 32)
		go _chunkify(path, readsPerChunk, ref, bai, chunks[ref.Name()])
	}

	merged := make(chan string)
	go func() {
		for {
			anyFound := false
			for chrom := range chunks {
				if chunks[chrom] == nil {
					continue
				}
				select {
				case chunk, ok := <-chunks[chrom]:
					if !ok {
						chunks[chrom] = nil
						continue
					}
					merged <- chunk
					anyFound = true
				case <-time.After(time.Millisecond):
					anyFound = true
				}
			}
			if !anyFound {
				close(merged)
				break
			}

		}
	}()
	return merged
}
func _chunkify(path string, readsPerChunk int, ref *sam.Reference, bai *bam.Index, ch chan string) {
	f, err := os.Open(path)
	if err != nil {
		close(ch)
		log.Fatal(err)
	}
	b, err := bam.NewReader(f, 0)
	if err != nil {
		close(ch)
		log.Fatal(err)
	}
	stats, ok := bai.ReferenceStats(ref.ID())
	if !ok {
		log.Printf("cant get stats for reference: %q. Assuming it has no reads\n", ref)
		close(ch)
		return
	}
	nReads := stats.Mapped
	if nReads == 0 {
		log.Printf("no reads mapped for: %q\n.", ref)
		close(ch)
		return
	}
	if nReads < uint64(1.2*float64(readsPerChunk)) {
		ch <- fmt.Sprintf("%s\t0\t%d\t%d", ref.Name(), ref.Len(), nReads)
		close(ch)
		return
	}

	log.Printf("chunking: %s. %d reads into %d / chunk \n", ref, nReads, readsPerChunk)

	it, err := bam.NewIterator(b, bai.Chunks(ref, 0, ref.Len()))
	if err != nil {
		log.Fatal(err)
	}

	cnt := 0
	chunk := 0
	lastStart := 0
	var aln *sam.Record

	for it.Next() {
		aln = it.Record()
		if aln.Flags&sam.Unmapped != 0 {
			continue
		}
		cnt++
		if cnt == readsPerChunk {
			ch <- fmt.Sprintf("%s\t%d\t%d\t\t%d", ref.Name(), lastStart, aln.Start(), cnt)
			cnt = 0
			lastStart = aln.Start()
			chunk++
		}
	}
	// leftover reads from reads modulo readsPerChunk
	if cnt > 0 {
		ch <- fmt.Sprintf("%s\t%d\t%d\t\t%d", ref.Name(), lastStart, aln.Start(), cnt)
	}
	err = it.Close()
	if err != nil {
		log.Println(err)
	}
	close(ch)
}

func main() {
	chunks := flag.Int("chunks", 1000, "approximate number of chunks to split bam file")
	flag.Parse()
	bam := flag.Arg(0)
	for chunk := range chunkify(bam, *chunks) {
		fmt.Println(chunk)
	}
}

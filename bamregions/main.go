package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"

	"github.com/brentp/xopen"
)

func main() {
	f := flag.Int("f", -1, "required sam flags")
	flag.Parse()
	if len(flag.Args()) < 2 {
		fmt.Fprintf(os.Stderr, "Usage of %s\n:", os.Args[0])
		fmt.Fprintf(os.Stderr, "\t%s {bed} {bam}", os.Args[0])
		flag.PrintDefaults()
		os.Exit(1)
	}

	bed := flag.Arg(0)
	bam := flag.Arg(1)
	n := runtime.GOMAXPROCS(-1)

	work := make(chan chan RdrBuf, n-1)
	go fillWork(work, bed, 100, bam, *f)

	wtr := exec.Command("sambamba", "view", "-f", "bam", "--sam-input", "/dev/stdin")
	wbam, err := wtr.StdinPipe()
	check(err)
	wtr.Stdout = os.Stdout
	wtr.Stderr = os.Stderr
	wtr.Start()

	for ch := range work {
		rdr := <-ch
		wbam.Write(rdr.buf)
		io.Copy(wbam, rdr.rdr)
		rdr.rdr.Close()
	}
	wbam.Close()

}

type RdrBuf struct {
	rdr io.ReadCloser
	buf []byte
}

func fillWork(work chan chan RdrBuf, bed string, nRegions int, bam string, flag int) {

	j := 0
	for regs := range Regions(bed, nRegions, 40) {
		args := []string{"view", bam}
		if j == 0 {
			args = append(args, "-h")
		}
		if flag != -1 {
			args = append(args, []string{"-f", strconv.Itoa(flag)}...)
		}
		j += 1
		for _, r := range regs {
			args = append(args, r.String())
		}
		ch := make(chan RdrBuf)
		work <- ch
		go func() {
			rdr, buf := GetUrl(bam, args)
			ch <- RdrBuf{rdr, buf}
			close(ch)
		}()

	}
	close(work)

}

type Region struct {
	chrom string
	start int
	end   int
}

func (r Region) String() string {
	return fmt.Sprintf("%s:%d-%d", r.chrom, r.start, r.end)
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func Regions(bed string, nRegions int, expand int) chan []Region {
	ch := make(chan []Region)

	go func() {
		rdr, err := xopen.Ropen(bed)
		check(err)
		defer rdr.Close()
		regions := make([]Region, 0, nRegions)
		for {
			line, err := rdr.ReadString('\n')
			if err == io.EOF {
				break
			} else {
				check(err)
			}
			toks := strings.Split(strings.TrimSuffix(line, "\n"), "\t")
			if toks[0][0] == '#' {
				continue
			}
			s, err := strconv.Atoi(toks[1])
			if s > expand {
				s -= expand
			}
			check(err)
			e, err := strconv.Atoi(toks[2])
			e += expand
			check(err)
			if s == e {
				continue
			}
			r := Region{toks[0], s, e}
			regions = append(regions, r)
			if len(regions) == nRegions {
				ch <- regions
				regions = make([]Region, 0, nRegions)
			}

		}
		if len(regions) > 0 {
			ch <- regions
		}
		close(ch)

	}()
	return ch
}

func GetUrl(url string, args []string) (io.ReadCloser, []byte) {
	args = append(args, url)
	cmd := exec.Command("samtools", args...)
	rpipe, err := cmd.StdoutPipe()
	check(err)
	cmd.Start()
	buf := make([]byte, 32*65536)
	n, err := rpipe.Read(buf)
	if err != io.EOF {
		check(err)
	}

	return rpipe, buf[:n]
}

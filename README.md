miscellaneous script-like stuff in go for bioinformatics

strip-sam-aux
-------------

remove select aux tags from a *bam* file (reads a bam, writes a stripped bam)

somatics
--------

call somatic variants on a multi-sample VCF with a single normal sample and many tumor samples (from the same patient, presumably across time).
A site where any tumor sample has a variant relative to the normal will have PASS or "." as the filter. This uses
the genotype-likelihood method from @chapmanb 's bcbio and @cc2qe 's speedseq, but handles multiple tumor samples.

chunkbam
--------

To parallelize genomic operations, we must split the genome. `chunkbam` splits the genome
as evenly as possible given an indexed bam and the number of desired chunks. It outputs
a BED file where **each region contains an equal number of mapped reads**. On a 7GB RNA-Seq
file (just as an example because it has very uneven coverage) it starts outputting rows
ready for parallelization after 6 seconds and finishes in under 1.5 minutes of user time
when using 20 CPUs.

Usage:

```Shell
	chunkbam -min-gap 10 -chunks 2000 some.bam > even-regions.bed
```

Then even-regions.bed will look like:
```
hs37d5	0	35477943	68644
Y	0	59373566	76386
13	0	27826971		108961
6	0	7570783		108961
20	0	2442466		108961
7	0	2674992		108961
10	0	5495221		108961
18	0	3458280		108961
MT	0	1686		108961
3	0	8607122		108961
```

Note that it attempts to get an even number (in this case *108961*) of reads per region,
but for smaller chromosomes it just takes whats available. There will also be fewer
reads in regions at the ends of chromsomes.

Also note that order is preserved only within chromosomes, but not across chromosomes due
to internal parallelization.

In many cases, we want to split when there are gaps in coverage. the `-min-gap` parameter allows this.
In the example above, once the chunksize is reached, a gap of 10 bases without coverage is also required
in order to make a new split.

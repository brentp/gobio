miscellaneous script-like stuff in go for bioinformatics

strip-sam-aux
-------------

remove select aux tags from a *bam* file (reads a bam, writes a stripped bam)

somatics
--------

call somatic variants on a multi-sample VCF with a single normal sample and many tumor samples (from the same patient, presumably across time).
A site where any tumor sample has a variant relative to the normal will have PASS or "." as the filter. This uses
the genotype-likelihood method from @chapmanb 's bcbio and @cc2qe 's speedseq, but handles multiple tumor samples.



# deplete-R

This is a module for scanning FASTQ files for sequencing reads which match a pattern, assembling a library of sequences and read counts for a degenerate loci next to the matched pattern, and then building position-weight matrixes based on the depletion or enrichment of those degenerate sequences in control vs experimental data sets.

# Objective

So, I've had occasion to use sequencing on the output of experiments like SELEX or plasmid depletion. In these experiments, degenerate loci are introduced into DNA constructs with degnerate oligos (usually), enriched or depleted in some way, and then the sequences which favor enrichment or depletion are resolved by short-read NGS.

Now, there are really strong workflows already built for the data analysis of these types of experiments. Sometimes though, I have a ton of different conditions and a huge pile of fastq files, and I know only a small proportion of them will have actually worked, the rest being failures for whatever reason (in a big screening experiment, they can't all be winners, right?)

So, this script was written to blitz through fastQ files, find the degenerate sites, and do quick-and-dirty PWMs for enrichment/depletion that can be easily viewed as Sequence Logos, and I can concentrate on the experiments that seem like they've produced a non-negative outcome. Are sequence logos perfect for this? Definately not, but that's why it's called quick and dirty, and not slow and perfect.

# Usage

It's pretty straightforward, give the pwmFromFastq function a vector of negative control files, a vector of experimental files, and the invariant sequence. The fastQ files do need to be deindexed prior to this though. It returns a list of two PWMs corresponding to enriched and depleted, and SeqLogo will take a PWM and make a logo with no additional information needed.

# Requirements

The ShortRead and seqLogo packages out of the Bioconductor library (and all their dependencies). Built with ShortRead v1.28.0 and seqLogo v1.36.0 using R v3.2 and RStudio v1.0.143

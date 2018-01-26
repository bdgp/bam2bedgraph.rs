cassette_reannotation
=====================

A set of rust tools used in the Conboy paper analysis.

Requires:

- latest Rust stable (or nightly)

To build:

Run ```cargo build```

Here is the list of tools and their description:

- bam2bedgraph

`Convert a BAM annotation file into a bedGraph or bigWig file.`

- cassette_reannotation

`Given an input annotation and a set of BAM files, discover any unannotated
cassette exons and produce a reannotation with new transcript features
that include the discovered cassettes.`

- cassette_lengths

`Produce a histogram of cassette exon lengths within the intron of a
transcript in another gene.`

- exon_cov

`Given an input annotation and a set of BAM files, produce a list of exon
base coverages and RPKM values.`

- adjusted_intron_psi

`Given an input annotation, a set of BAM files, and a SplAdder retained intron
output table, produce a new retained intron table with corrected
adjusted PSI values.`

More information about each tool can be found by executing the tool
with the `--help` argument.
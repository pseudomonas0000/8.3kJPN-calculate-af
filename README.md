# 8.3kJPN-calclate-af 
Coded by pseudomonas0000, 2021/3/18
Tested on macOS Catalina (10.15.7)

## Calculate allele frequency from 8.3kJPN Genotype datasets to enable anotation
The ToMMo [8.3KJPN Genotype Frequency Panel](https://jmorp.megabank.tohoku.ac.jp/202102/downloads/#variant) datasets have some INFO tags. This inclueds allele frequency and genotype information. In oder to obrain accurate allele frequency, it must be recalculated from genotype information. This script provides accurate allele freqyency and be using for annotation your vcf file.

## Requirements
* user/bin/bash
* Perl (user/bin/perl)
* Some tools (bcftools, bgzip, tabix, vt)

## Procedure
1. Download this repository and the [8.3kJPN Genotype datasets](https://jmorp.megabank.tohoku.ac.jp/202102/downloads/#variant) (four `.gz` files and four `.tbi` files). Save the file in the repository.
1. On your console, run `bash build.sh path/to/hg19.fasta output-file-name` . The bash script has two command line aurguments, first is full path to refseq hg19 fasta file, second is an output file name.
1. Annotate the created file after runnning. On your console, `snpsift annotate created-file.vcf your.vcf > annotated.vcf` .

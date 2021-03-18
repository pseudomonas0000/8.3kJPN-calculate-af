#!/usr/bin/bash

# download datasets from jmorp website site (login needed and wget command can not use.)
	# tommo-8.3kjpn-20200831-gf_snvall-autosome.vcf.gz
	# tommo-8.3kjpn-20200831-gf_indelall-autosome.vcf.gz	tommo-8.3kjpn-20200831-gf_snvall-autosome.vcf.gz.tbi
	# tommo-8.3kjpn-20200831-gf_indelall-autosome.vcf.gz.tbi	tommo-8.3kjpn-20200831-gf_snvall-chrX_PAR3.vcf.gz
	# tommo-8.3kjpn-20200831-gf_snvall-chrX_PAR3.vcf.gz.tbi	tommo-8.3kjpn-20200831-gf_indelall-chrX_PAR3.vcf.gz.tbi


# concatate 4 files (autosomal snv, autosomal indel, chrX snv and chrX indel) and indexing
bcftools concat -Ou -a tommo-8.3kjpn-20200831-gf_snvall-autosome.vcf.gz tommo-8.3kjpn-20200831-gf_indelall-autosome.vcf.gz tommo-8.3kjpn-20200831-gf_indelall-chrX_PAR3.vcf.gz tommo-8.3kjpn-20200831-gf_snvall-chrX_PAR3.vcf.gz | bcftools sort - -O z > tommo-8.3kjpn-all.vcf.gz
tabix -s 1 -b 2 -e 2 tommo-8.3kjpn-all.vcf.gz

# transfer annotation from INFO column to FORMAT.
	# For convenience, add the sample name (in this case  "tommo") to the header.
	bcftools view -h tommo-8.3kjpn-all.vcf.gz |gsed "s/FILTER\tINFO/FILTER\tINFO\tFORMAT\ttommo/" > new_header.txt
	bcftools reheader -h new_header.txt tommo-8.3kjpn-all.vcf.gz > tommo-8.3kjpn-all_reheader.vcf.gz
	
	# Extract genotype information in INFO column to a tab-delimited annotation file. 
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%GenotypeCount\t%GenotypeCount_MALE\t%GenotypeCount_FEMALE\n' tommo-8.3kjpn-all_reheader.vcf.gz | bgzip -c > annot.txt.gz
	# Index the file with tabix.
	tabix -s1 -b2 -e2 annot.txt.gz
	
	# Add a header information for new annotation, transfer the annotation to sample 'tommo', add CHROM column to 'chr'and remove unnecessary INFO tags.
	bcftools annotate --rename-chrs jmorp-annotation/rename-chr.txt -s tommo -a annot.txt.gz -h jmorp-annotation/header-add1.txt -c CHROM,POS,REF,ALT,FMT/GenotypeCount,FMT/GenotypeCount_MALE,FMT/GenotypeCount_FEMALE -x ^INF/GenotypeCount,INF/GenotypeCount_MALE,INF/GenotypeCount_FEMALE -O v tommo-8.3kjpn-all_reheader.vcf.gz > tommo-8.3kjpn-all.genotype.vcf

# Decompose and normalize ALT column.
vt decompose -s tommo-8.3kjpn-all.genotype.vcf | vt normalize - -r ~/Downloads/hg19/hg19.fasta -o tommo-8.3kjpn-all.decomp.vcf

# Calclate allele frequenies from genotype counts, etc.
perl jmorp-annotation/calc-af.pl tommo-8.3kjpn-all.decomp.vcf
echo "calc-af.pl finished"

# For calclated allele frequencies, add header lines, then remove unnecesary format column and info tags.
bcftools annotate -h jmorp-annotation/header-add2.txt -x FMT,INF/OLD_MULTIALLELIC,INF/OLD_VARIANT,INF/GenotypeCount,INF/GenotypeCount_MALE,INF/GenotypeCount_FEMALE tommo-8.3kjpn-all.decomp.calc.vcf  -O z -o tommo-8.3kjpn-all.decomp.calc.rminfo.vcf.gz

echo "Building completed"
date

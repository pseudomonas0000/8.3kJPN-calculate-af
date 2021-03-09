# download datasets from jmorp site
	# tommo-8.3kjpn-20200831-gf_snvall-autosome.vcf.gz
	# tommo-8.3kjpn-20200831-gf_indelall-autosome.vcf.gz	tommo-8.3kjpn-20200831-gf_snvall-autosome.vcf.gz.tbi
	# tommo-8.3kjpn-20200831-gf_indelall-autosome.vcf.gz.tbi	tommo-8.3kjpn-20200831-gf_snvall-chrX_PAR3.vcf.gz
	# tommo-8.3kjpn-20200831-gf_snvall-chrX_PAR3.vcf.gz.tbi	tommo-8.3kjpn-20200831-gf_indelall-chrX_PAR3.vcf.gz.tbi

# 3つのvcfファイルを結合しソートする autosomal indel, chrX snv, ChrX indel
# 常染色体snvはサイズが大きいため、のちに個別に処理する -Ouはパイプでbcftoolsのサブコマンドをつなぐ時につける 高速になる
bcftools concat  -Ou -a tommo-8.3kjpn-20200831-gf_indelall-autosome.vcf.gz tommo-8.3kjpn-20200831-gf_indelall-chrX_PAR3.vcf.gz tommo-8.3kjpn-20200831-gf_snvall-chrX_PAR3.vcf.gz | bcftools sort - -O z > tommo-8.3kjpn-concat.sort.vcf.gz
	bcftools concat -a tommo-8.3kjpn-20200831-gf_indelall-autosome.vcf.gz tommo-8.3kjpn-20200831-gf_indelall-chrX_PAR3.vcf.gz tommo-8.3kjpn-20200831-gf_snvall-chrX_PAR3.vcf.gz -o tommo-8.3kjpn-concat.vcf.gz

	# 本番用
	# concatファイルに対してindexファイルを作る前にsortをする -Oは出力ファイルのタイプ指定 zは圧縮vcfファイル
	bcftools sort tommo-8.3kjpn-concat.vcf.gz -O z >tommo-8.3kjpn-concat.sort.vcf.gz
# 必要のないINFO情報は消し、CHROM列に「chr」をつける 
bcftools annotate -Ou -x ^INF/GenotypeCount,INF/GenotypeCount_MALE,INF/GenotypeCount_FEMALE tommo-8.3kjpn-concat.sort.vcf.gz | bcftools annotate --rename-chrs rename-chr.txt - -O z -o tommo-8.3kjpn-concat.sort.rminfo-chr.vcf.gz

# 常染色体snvに対してもINFO情報の削除とchrをつける
bcftools annotate -Ou -x ^INF/GenotypeCount,INF/GenotypeCount_MALE,INF/GenotypeCount_FEMALE tommo-8.3kjpn-20200831-gf_snvall-autosome.vcf.gz | bcftools annotate --rename-chrs rename-chr.txt - -O z -o tommo-8.3kjpn-snvall-autosome.rminfo-chr.vcf.gz

# concatするのでindexをつける (複数ファイルまとめてできそう)
	bcftools index tommo-8.3kjpn-concat.sort.rminfo-chr.vcf.gz
	bcftools index tommo-8.3kjpn-snvall-autosome.rminfo-chr.vcf.gz

# 常染色体snvファイルとその他のファイルを結合してsortする -O vは非圧縮vcf形式で出力 數十分かかった
bcftools concat -Ou -a tommo-8.3kjpn-snvall-autosome.rminfo-chr.vcf.gz tommo-8.3kjpn-concat.sort.rminfo-chr.vcf.gz | bcftools sort - -O v > tommo-8.3kjpn-full.vcf

# perlを実行 ファイル名に注意
perl move-genotype.pl tommo-8.3kjpn-full.vcf
	# headerの最後に「FORMAT」とサンプル名を追加 (サンプル名の付加にはbcftools reheader -s, --samplesを試したがダメだった)
	# autosomeとchrXそれぞれにINFO列内のGenotype情報をGenotype列に移動した

	# bcftoolsでもできそう samtools.github.io/bcftools/howto
	# Extract INFO/GenotypeCount, INFO/GenotypeCount_MALE, INFO/GenotypeCount_FEMALE into a tab-delimited annotation file
	# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/GenotypeCount'

# vcfファイルのヘッダーを追加する (header-add1.txt)　数十分かかった
bcftools annotate -h header-add1.txt tommo-8.3kjpn-full.move.vcf >tommo-8.3kjpn-full.move.header1.vcf

# (必要なら)decompose前にソートなどの確認
# snpsift vcfcheck

# vtでdecomposeとnormalize 数十分かかった
vt decompose -s *header1.vcf | vt normalize - -r ~/Desktop/panel_bed/hg19.fasta -o tommo-8.3kjpn-full.move.header1.decomp-normal.vcf
	# -sをつけないとGenotype fieldsのGT,PL,GL,DP以外の情報が失われる。type=Gだけではなく、type=Aも失われる。

	# decompose前後
	# chr1	103379828	rs145342760	AT	ATT,A	.	PASS	AC=342,5;AN=16760;AF=0.0204,0.0003;AC_MALE=166,0;AN_MALE=7834;AF_MALE=0.0212,0;AC_FEMALE=176,5;AN_FEMALE=8926;AF_FEMALE=0.0197,0.0006;TOMMO_PLATFORM_BIAS_TEST_PVALUE=0.1474	GenotypeCount:GenotypeCount_MALE:GenotypeCount_FEMALE	8037,334,4,5,0,0:3754,160,3,0,0,0:4283,174,1,5,0,0

	# chr1	103379828	rs145342760	AT	ATT	.	PASS	AC=342;AN=16760;AF=0.0204;AC_MALE=166;AN_MALE=7834;AF_MALE=0.0212;AC_FEMALE=176;AN_FEMALE=8926;AF_FEMALE=0.0197;TOMMO_PLATFORM_BIAS_TEST_PVALUE=0.1474;OLD_MULTIALLELIC=chr1:103379828:AT/ATT/A	GenotypeCount:GenotypeCount_MALE:GenotypeCount_FEMALE	8037,334,4:3754,160,3:4283,174,1
	# chr1	103379828	rs145342760	AT	A	.	PASS	AC=5;AN=16760;AF=0.0003;AC_MALE=0;AN_MALE=7834;AF_MALE=0;AC_FEMALE=5;AN_FEMALE=8926;AF_FEMALE=0.0006;TOMMO_PLATFORM_BIAS_TEST_PVALUE=0.1474;OLD_MULTIALLELIC=chr1:103379828:AT/ATT/A	GenotypeCount:GenotypeCount_MALE:GenotypeCount_FEMALE	8037,5,0:3754,0,0:4283,5,0

# perlスクリプトでgenotype列から計算を行い、INFO列に入れる (数時間?かかった)
perl calc-af.pl *decomp-normal.vcf

# headerを追加する(header-add2.txt)、必要のないINFO情報とFORMAT情報を除く
bcftools annotate -h header-add2.txt -x FMT,INF/OLD_MULTIALLELIC,INF/OLD_VARIANT,INF/GenotypeCount,INF/GenotypeCount_MALE,INF/GenotypeCount_FEMALE *.calc.vcf -o tommo-8.3kjpn-full.move.header1.decomp-normal.calc.rminfo.vcf

# パネル領域のみを取り出す
bedtools intersect -header -a tommo-8.3kjpn-full.move.header1.decomp-normal.calc.rminfo.vcf -b all_panel.bed > tommo-8.3kjpn-full.move.header1.decomp-normal.calc.rminfo.allpanel210205.vcf

# 山口さんのファイルと見比べるため、必要な情報のみtabテキスト形式で出力
	bgzip -c tommo-8.3kjpn-full.move.header1.decomp-normal.calc.rminfo.allpanel210205.vcf > tommo-8.3kjpn-full.move.header1.decomp-normal.calc.rminfo.allpanel210205.vcf.gz
	bcftools index tommo-8.3kjpn-full.move.header1.decomp-normal.calc.rminfo.allpanel210205.vcf.gz
	bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/tommo8.3k_AC\t%INFO/tommo8.3k_AF\t%INFO/tommo8.3k_AN\t%INFO/tommo8.3k_exp_hom\t%INFO/tommo8.3k_carrier_het\t%INFO/tommo8.3k_hom\t%INFO/tommo8.3k_carrier_hom\t%INFO/tommo8.3k_hemi\t%INFO/tommo8.3k_carrier_hemi\n' tommo-8.3kjpn-full.move.header1.decomp-normal.calc.rminfo.allpanel210205.vcf.gz -H -o test.txt



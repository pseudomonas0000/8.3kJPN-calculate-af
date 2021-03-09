#!/usr/bin/perl

use strict;
use warnings;

# 使用するモジュールと関数
use File::Copy qw/copy move/;

my $infile = $ARGV[0];
my $tmpfile = "tmp$$";

$infile =~ /(.+)\.vcf/;
my $outfile = $1 . ".calc.vcf";

open my $infh, "<", $infile
	or die qq/Cannot open $infile : $!/;
open my $tmpfh, ">", $tmpfile
	or die qq/Cannot open $tmpfile : $!/;

sub auto_calc{
# [all_rr, all_ra, all_aa]
	my @genotype_auto = @_;

	my $auto_ac = (2 * $genotype_auto[2]) + $genotype_auto[1];
	my $auto_an = 2 * ($genotype_auto[0] + $genotype_auto[1] + $genotype_auto[2]);
	
	# calculate carrier hom
	undef(my $auto_carrier_hom);
	if($genotype_auto[2] != 0){
		$auto_carrier_hom = int(($auto_an / 2) / $genotype_auto[2]);
	}else{
		$auto_carrier_hom = ".";
	}

	# calculation carrier het
	undef(my $auto_carrier_het);
	if($genotype_auto[1] != 0){
		$auto_carrier_het = int(($auto_an / 2) / $genotype_auto[1]);
	}else{
		$auto_carrier_het = ".";
	}
	
	# calcutation af
	undef(my $auto_af);
	if($auto_an != 0){
		$auto_af = sprintf("%.4g", ($auto_ac / $auto_an));
	}else{
		$auto_af = ".";
	}

	# calculate expected hom
	undef(my $auto_exp_hom);
	if($genotype_auto[1] != 0){
		$auto_exp_hom = ((($auto_an / 2) / $genotype_auto[1]) ** 2) * 4;
	}else{
		$auto_exp_hom = ".";
	}

	# @array_info["tommo8.3k_AC", "tommo8.3k_AN", "tommo8.3k_AF", "tommo8.3k_carrier_het", "tommo8.3k_hom", "tommo8.3k_carrier_hom", "tommo8.3k_exp_hom", "tommo8.3k_hemi", "tommo8.3k_carrier_hemi"]
	my $result_auto_calc =[$auto_ac, $auto_an, $auto_af, $auto_carrier_het, $genotype_auto[2], $auto_carrier_hom, $auto_exp_hom, ".", "."];
	return $result_auto_calc;
}

sub chrx_calc1{
# [$male_r, male_a, female_rr, female_ra, female_aa]
	my @count = @_;

	my $all_ac = $count[1] + $count[3] + (2 * $count[4]);
	my $male_an = $count[0] + $count[1];
	my $female_an = 2 * ($count[2] + $count[3] + $count[4]);
	my $all_an = $male_an + $female_an;
	my $male_hemi = $count[1];
	my $female_hom = $count[4];

	# calculation all_af
	my $all_af;
	if($all_an != 0){
		$all_af = sprintf("%.4g", ($all_ac / $all_an));
	}else{
		$all_af = ".";
	}

	# calculation $male_carrier_hemi
	my $male_carrier_hemi;
	if($male_hemi != 0){
		$male_carrier_hemi = int($male_an / $male_hemi);
	}else{
		$male_carrier_hemi = ".";
	}

	# calculation $female_carrier_hom
	my $female_carrier_hom;
	if($female_hom != 0){
		$female_carrier_hom = int(($female_an / 2) / $female_hom);
	}else{
		$female_carrier_hom = ".";
	}

	# @array_info["tommo8.3k_AC", "tommo8.3k_AN", "tommo8.3k_AF", "tommo8.3k_carrier_het", "tommo8.3k_hom", "tommo8.3k_carrier_hom", "tommo8.3k_exp_hom", "tommo8.3k_hemi", "tommo8.3k_carrier_hemi"]
	my $result_chrx_calc = [$all_ac, $all_an, $all_af, ".", $female_hom, $female_carrier_hom, ".", $male_hemi, $male_carrier_hemi];
	return $result_chrx_calc;
}

sub chrx_calc2{
# pseudoautosomal regions (PAR)
# [male_rr, male_ra, male_aa, female_rr, female_ra, female_aa]
	my @count = @_;

	my $all_ac = $count[1] + $count[4] + 2 * ($count[2] + $count[5]);
	my $all_an = 2 * ($count[0] + $count[1] + $count[2] + $count[3] + $count[4] + $count[5]);

	# calculation all_af
	my $all_af;
	if($all_an != 0){
		$all_af = sprintf("%.4g", ($all_ac / $all_an));
	}else{
		$all_af = ".";
	}

	# calculation carrier_het (male + female)
	my $all_het = $count[1] + $count[4];
	my $all_carrier_het;
	if( $all_het != 0){
		$all_carrier_het = int(($all_an / 2) / $all_het);
	}else{
		$all_carrier_het = ".";
	}

	# calculation carrier_hom (male + female)
	my $all_hom = $count[2] + $count[5];
	my $all_carrier_hom;	
	if($all_hom != 0){
		$all_carrier_hom = int(($all_an / 2) / $all_hom);
	}else{
		$all_carrier_hom = ".";
	}

	# @array_info["tommo8.3k_AC", "tommo8.3k_AN", "tommo8.3k_AF", "tommo8.3k_carrier_het", "tommo8.3k_hom", "tommo8.3k_carrier_hom", "tommo8.3k_exp_hom", "tommo8.3k_hemi", "tommo8.3k_carrier_hemi"]
	my $result_chrx_calc2 = [$all_ac, $all_an, $all_af, $all_carrier_het, $count[5], $all_carrier_hom, ".", ".", "."];
	return $result_chrx_calc2;
}

while(<$infh>){
	my $inline = $_;
	chomp($inline);
	my @inarray = split(/\t/, $inline);

	if($inarray[0] =~ /^#/){
		if($inarray[0] =~ /^#CHROM/){
			# CHROM行のサンプル名とFORMATを消して出力
			print $tmpfh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
			next;
		}else{
			print $tmpfh "$inline\n";
			next;	
		}
	}

	undef(my @array_calc);
	undef(my @array_info);
	@array_info = ("tommo8.3k_AC", "tommo8.3k_AN", "tommo8.3k_AF", "tommo8.3k_carrier_het", "tommo8.3k_hom", "tommo8.3k_carrier_hom", "tommo8.3k_exp_hom", "tommo8.3k_hemi", "tommo8.3k_carrier_hemi");
	
	# autosome
	if($inarray[0] =~ /chr\d+/){

		undef(my $auto_calc_ref);

		$inarray[9] =~ /(\d+)\,(\d+)\,(\d+):(\d+)\,(\d+)\,(\d+):(\d+)\,(\d+)\,(\d+)/;
		$auto_calc_ref = auto_calc($1, $2, $3, $4, $5, $6, $7, $8, $9);
		@array_calc = @{$auto_calc_ref};

	# chrX
	}elsif($inarray[0] eq "chrX"){

		undef(my $chrx_calc_ref);
		
		# pseudoautosomal region (PAR) or not
		if($inarray[9] =~ /^(\d+)\,(\d+):(\d+)\,(\d+)\,(\d+)/){			
			$chrx_calc_ref = chrx_calc1($1, $2, $3, $4, $5);
		}elsif($inarray[9] =~ /^(\d+)\,(\d+)\,(\d+):(\d+)\,(\d+)\,(\d+)/){
			$chrx_calc_ref = chrx_calc2($1, $2, $3, $4, $5, $6);
		}
		@array_calc = @{$chrx_calc_ref};
	}

	for (my $i = 0; $i < @array_calc; $i++){
		$array_info[$i] .= "=" . $array_calc[$i];
	}
	my $line_info = ";" . join(";", @array_info);
	$inarray[7] .= $line_info;

	my $outline = join("\t", @inarray);
	print $tmpfh "$outline\n";
}

close($infh);
close($tmpfh);

move $tmpfile, $outfile 
	or die qq/Cannot move $tmpfh to $outfile :$1/;


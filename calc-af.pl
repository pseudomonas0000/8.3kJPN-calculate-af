#!/usr/bin/perl

use strict;
use warnings;

# use module and methods
use File::Copy qw/copy move/;
use List::Util qw/sum/;

my $infile = $ARGV[0];
my $tmpfile = "tmp$$";
(my $outfile = $infile) =~ s/.vcf/.calc.vcf/;

open my $infh, "<", $infile
	or die qq/Cannot open $infile : $!/;
open my $tmpfh, ">", $tmpfile
	or die qq/Cannot open $tmpfile : $!/;

# autosome and pseudoautosomal regions
sub auto_calc{
# the contents of the @_ are [male_rr, male_ra, male_aa, female_rr, female_ra, female_aa]
	my @auto_count = @_;

	my $geno_hom = $auto_count[2] + $auto_count[5];
	my $geno_het = $auto_count[1] + $auto_count[4];
	my $auto_ac = (2 * $geno_hom) + $geno_het;
	my $auto_an = 2 * sum(@auto_count);
	
	# calculate carrier hom
	undef(my $auto_carrier_hom);
	if($geno_hom != 0){
		$auto_carrier_hom = int(($auto_an / 2) / $geno_hom);
	}else{
		$auto_carrier_hom = ".";
	}

	# calculate carrier het
	# calculate expected hom	
	undef(my $auto_carrier_het);
	undef(my $auto_exp_hom);
	if($geno_het != 0){
		$auto_carrier_het = int(($auto_an / 2) / $geno_het);
		$auto_exp_hom = int(((($auto_an / 2) / $geno_het) ** 2) * 4);
	}else{
		$auto_carrier_het = ".";
		$auto_exp_hom = ".";
	}
	
	# calcutation af
	undef(my $auto_af);
	if($auto_an != 0){
		$auto_af = sprintf("%.4f", ($auto_ac / $auto_an));
	}else{
		$auto_af = ".";
	}

	# @array_info = ("tommo8.3k_AC", "tommo8.3k_AN", "tommo8.3k_AF", "tommo8.3k_carrier_het", "tommo8.3k_hom", "tommo8.3k_carrier_hom", "tommo8.3k_exp_hom", "tommo8.3k_hemi", "tommo8.3k_carrier_hemi");
	my $result_auto_calc = [$auto_ac, $auto_an, $auto_af, $auto_carrier_het, $geno_hom, $auto_carrier_hom, $auto_exp_hom, ".", "."];
	return $result_auto_calc;
}

# chromosome X
sub chrx_calc{
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
		$all_af = sprintf("%.4f", ($all_ac / $all_an));
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

	# @array_info = ("tommo8.3k_AC", "tommo8.3k_AN", "tommo8.3k_AF", "tommo8.3k_carrier_het", "tommo8.3k_hom", "tommo8.3k_carrier_hom", "tommo8.3k_exp_hom", "tommo8.3k_hemi", "tommo8.3k_carrier_hemi");
	my $result_chrx_calc = [$all_ac, $all_an, $all_af, ".", $female_hom, $female_carrier_hom, ".", $male_hemi, $male_carrier_hemi];
	return $result_chrx_calc;
}

my $j = 0;
while(<$infh>){
	my $inline = $_;
	chomp($inline);
	my @inarray = split(/\t/, $inline);

	if($inarray[0] =~ /^#/){
		if($inarray[0] =~ /^#CHROM/){
			# remove 'FORMAT' and 'tommo'.
			print $tmpfh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
		}else{
			print $tmpfh "$inline\n";
		}
		next;
	}

	undef(my @array_calc);
	undef(my @array_info);
	@array_info = ("tommo8.3k_AC", "tommo8.3k_AN", "tommo8.3k_AF", "tommo8.3k_carrier_het", "tommo8.3k_hom", "tommo8.3k_carrier_hom", "tommo8.3k_exp_hom", "tommo8.3k_hemi", "tommo8.3k_carrier_hemi");
	
	undef(my $calc_ref);
	# autosome and pseudoautosomal region (PAR)
	if($inarray[9] =~ /:?(\d+)\,(\d+)\,(\d+):(\d+)\,(\d+)\,(\d+)$/){
		$calc_ref = auto_calc($1, $2, $3, $4, $5, $6);
	# chrX 
	}elsif($inarray[9] =~ /^[\.:]*(\d+)\,(\d+):(\d+)\,(\d+)\,(\d+)/){
		$calc_ref = chrx_calc($1, $2, $3, $4, $5, $6);
	}
	
	@array_calc = @{$calc_ref};

	for(my $i = 0; $i < @array_calc; $i++){
		$array_info[$i] .= "=" . $array_calc[$i];
	}
	my $line_info = ";" . join(";", @array_info);
	$inarray[7] .= $line_info;

	my $outline = join("\t", @inarray);
	print $tmpfh "$outline\n";

	print STDERR "READ\t$j\t$inarray[0]\t$inarray[1]\n" if($j % 1000000 == 0);
	$j ++;
}

close($infh);
close($tmpfh);

move $tmpfile, $outfile 
	or die qq/Cannot move $tmpfh to $outfile :$1/;


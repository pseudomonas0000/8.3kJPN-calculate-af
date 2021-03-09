#!/usr/bin/perl

use strict;
use warnings;

# 使用するモジュールと関数
use File::Copy qw/copy move/;

my $infile = $ARGV[0];
my $tmpfile = "move_genotype_tmp$$";

$infile =~ /(.+)\.vcf/;
my $outfile = $1 . ".move.vcf";

open my $infh, "<", $infile
	or die qq/Cannot open $infile : $!/;
open my $tmpfh, ">", $tmpfile
	or die qq/Cannot open $tmpfile : $!/;

while(<$infh>){
	my $inline = $_;
	chomp($inline);
	my @inarray = split(/\t/, $inline);

	if($inarray[0] =~ /^#/){
		if($inarray[0] =~ /^#CHROM/){
			print $tmpfh "$inline\tFORMAT\ttommo-gc\n";
			next;
		}else{
			print $tmpfh "$inline\n";
			next;	
		}
	}
	if($inarray[0] =~ /(chr\d+)/){
		$inarray[7] =~ /GenotypeCount=([\d\,]+);GenotypeCount_MALE=([\d\,]+);GenotypeCount_FEMALE=([\d\,]+)/;	
		print $tmpfh "$inline\tGenotypeCount:GenotypeCount_MALE:GenotypeCount_FEMALE\t$1:$2:$3\n";
	}elsif($inarray[0] eq "chrX"){
		$inarray[7] =~ /GenotypeCount_MALE=([\d\,]+);GenotypeCount_FEMALE=([\d\,]+)/;
		print $tmpfh "$inline\tGenotypeCount_MALE:GenotypeCount_FEMALE\t$1:$2\n";	
	}
}

close($infh);
close($tmpfh);

move $tmpfile, $outfile 
	or die qq/Cannot move $tmpfh to $outfile :$1/;
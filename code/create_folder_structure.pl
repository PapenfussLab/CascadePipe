#!/usr/bin/perl

use strict;

open(D,"../list.patients") || die "$!\n";

system("mkdir VCF");
system("mkdir VCF/pmcc_pipeline VCF/mpileup VCF/multiSNV");
system("mkdir TSV");

while(<D>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	my $case = $line[0];
		
	system("mkdir VCF/pmcc_pipeline/$case VCF/mpileup/$case VCF/multiSNV/$case");
	system("mkdir TSV/$case");
}
close(D);

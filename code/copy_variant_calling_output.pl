#!/usr/bin/perl

use strict;
use Cwd;

open(D,"../list.patients") || die "$!\n";

my $parent_dir = getcwd();

sub move
{
	my $files_ref = shift;
	my $dir = shift;

	my @files = @$files_ref;
		
	for(my $i=0;$i<@files;$i++)
	{
		chomp($files[$i]);
		system("mv $files[$i] $dir");
		print "mv $files[$i] $dir\n";
	}
}

while(<D>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	my $case = $line[0];

	print $case,"\n";	
	my @mpileups = `ls $parent_dir/../mpileup/$case\*\.mpileup\*vcf`;
	move(\@mpileups,"$parent_dir/VCF/mpileup\/$case\/");

	my @multisnvs = `ls $parent_dir/../multisnv/$case\_chr\*\_multiSNV.vcf`;
	move(\@multisnvs,"$parent_dir/VCF/multiSNV\/$case\/");

	my @pmccs = `ls $parent_dir/../seqliner/$case\*/Vcf/\*_merged_comb.vcf`;
	move(\@pmccs,"$parent_dir/VCF/pmcc_pipeline\/$case\/");

	my @tsvs = `ls $parent_dir/../seqliner/$case\*/Variants/\*_merged_comb.tsv`;
	move(\@tsvs,"$parent_dir\/TSV\/$case\/");
}
close(D);

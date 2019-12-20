#!/usr/bin/perl

use strict;
use Cwd;

my @chroms = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

open(D,$ARGV[0]) || die "$!\n";

my %case2bam=();

my $path = getcwd();

my %seen=();

while(<D>)
{
	chomp($_);
	next if ($_=~/^\#/);
	my @line = split(/\t/,$_);
	$line[0]=~/^((CAS|CA00)\d\d)/;
	my $case = $1;

	if(!defined $seen{$line[1]})
	{
		$seen{$line[1]}++;
		$case2bam{$case}=$line[1] .' ' . $case2bam{$case};
	}
	
	$case2bam{$case}.=$line[2] . ' ';
}
close(D);

for my $key (sort keys %case2bam)
{
	my @samples = split(/\s+/,$case2bam{$key});
	
	my $n = scalar(@samples);
	
	for(my $i=0;$i<@chroms;$i++)
	{
		my $command;
		$command="multiSNV -N$n --bam $case2bam{$key} --fasta /data/reference/indexes/human/g1k_v37/fasta/human_g1k_v37.fasta ";
                $command.= "--regions $chroms[$i] -d 10 -f $key\_chr$chroms[$i]\_multiSNV.vcf";


		open(H,"header.sh") || die "$!\n";

		my $out_file = "$key\_chr$chroms[$i]\_multiSNV" . '_run_multiSNV.sh';

		open(O,">$out_file");

		while(<H>)
		{
			chomp($_);
			my $job_name = "$key\_chr$chroms[$i]\_multiSNV";
			$job_name=~s/\-/\_/g;
			$_=~s/TEST/$job_name/;
			print O $_,"\n";
		}
		close(H);

		print O $command,"\n";
		close(O);


	}
	print "\n";
}



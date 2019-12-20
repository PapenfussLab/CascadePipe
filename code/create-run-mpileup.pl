#!/usr/bin/perl

use strict;

open(D,$ARGV[0]) || die "$!\n";
#has pairs of normal/tumor bam files

my $bed_file=$ARGV[1];

while(<D>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	
	my $out_file = $line[0] . ".mpileup";

	my $command1 = "samtools mpileup -l $bed_file -d 1000000 -A -B -f /data/reference/indexes/human/g1k_v37/fasta/human_g1k_v37.fasta $line[1] $line[2] \> $out_file";

	my $command2 = "VarScan mpileup2snp $out_file  --output-vcf 1 --min-avg-qual 20 --min-coverage 10 --min-reads2 2  --min-var-freq 0 --p-value 1 --strand-filter 1 > $out_file\.snp\.vcf";

	my $command3 = "VarScan mpileup2indel $out_file  --output-vcf 1 --min-avg-qual 20 --min-coverage 10 --min-reads2 2  --min-var-freq 0 --p-value 1 --strand-filter 1 > $out_file\.indel\.vcf";

	my $gzipped = $out_file . '.gz';

	my $command4 = "gzip $out_file";

	open(H,"header.sh") || die "$!\n";

	my $job_name = $line[0] . '_mpileup';
	$job_name=~s/\-/\_/g;

	my $bash_file = $line[0] . '_run_mpileup.sh';

	open(O,">$bash_file");

	while(<H>)
	{
		chomp($_);
		$_=~s/TEST/$job_name/;
		print O $_,"\n";
	}
	close(H);
	
	print O $command1,"\n";
	print O $command2,"\n";
	print O $command3,"\n";
	print O $command4,"\n";
}
close(D);


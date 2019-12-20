#!/usr/bin/perl

use strict;
use POSIX qw(ceil floor);

#GT:AD:BQ:DP:FA:SS       0/1:18,15:36:33:0.455:2


open(MAT,$ARGV[0]) || die "$!\n";

my $case = $ARGV[1];

my @impact_files = `ls /papenfuss_lab/analysis/CASCADE/FINAL/SNVs/TSV/$case/\*.tsv `;

my %mutation2gmaf=();
my %mutation2evsaafall=();


for(my $i=0;$i<@impact_files;$i++)
{
        chomp($impact_files[$i]);
        open(D,$impact_files[$i]) || die "$impact_files[$i]\: $!\n";

	my @name_info = split(/\//,$impact_files[$i]);
	$name_info[$#name_info]=~/(\S+)\_merged/;

	my $sample = $1;

        print STDERR "Processing file $impact_files[$i]\n";
	print STDERR "$sample\n";
	
	my $index_gmaf=-1;
	my $index_evsaafall=-1; 

        while(<D>)
        {
		my @line = split(/\t/,$_);
		if($_=~/^\#CHROM/)
		{
			for(my $i=0;$i<@line;$i++)
			{
                                if($line[$i] eq 'GMAF')
                                {
                                        $index_gmaf = $i;
                                }                       
                                if($line[$i] eq 'evs_aaf_all')
                                {
                                        $index_evsaafall=$i;
                                }
			}
			next;
		}
                next if ($_=~/^\#/);

		last if ($index_gmaf == -1 | $index_evsaafall == -1);

                my $key = "$line[0]\t$line[1]\t$line[3]\t$line[4]";

		$mutation2gmaf{$key} = $line[$index_gmaf];
		$mutation2evsaafall{$key} = $line[$index_evsaafall];
		
        }
        close(D);
}


while(<MAT>)
{
	chomp($_);
	if($_=~/sample/)
	{
		print $_,"\n";
		next;
	}
	else
	{
		my @line = split(/\t/,$_);
		my $key = "$line[2]\t$line[3]\t$line[4]\t$line[5]";

                if ($mutation2evsaafall{$key}=~/\d+/ || $mutation2gmaf{$key}=~/\d+/)
                {
                        print STDERR $key," excluded with GMAF ", $mutation2gmaf{$key}," and EVS ", $mutation2evsaafall{$key},"\n";
                        next;
                }
                else
                {
                       # print STDERR "KEPT $key\t$mutation2evsaafall{$key}\t$mutation2gmaf{$key}\n";
                        print $_,"\n";
                }
	}
}
close(MAT);

#!/usr/bin/perl

use strict;

#perl 3_create-dummy-vcf-for-multiSNV.pl ./VCF/multiSNV/CA-53/CA-53_multiSNV_matrix.txt


open(D,"../list.patients") || die "$!\n";

while(<D>)
{
        chomp($_);
        my @line = split(/\t/,$_);
        
        my $case = $line[0];
        print STDERR $case,"\n";

	my $matrix = "./VCF/multiSNV/$case/$case\_multiSNV_matrix.txt";

	open(MERGED_MAT,$matrix) || die "$!\n";
	my @names;

	system("rm ./VCF/multiSNV/$case\/\*_multiSNV_sample.vcf");

	while(<MERGED_MAT>)
	{
		chomp($_);
		if($_=~/var\tflag\t(.+)/)
		{
			@names = split(/\t/,$1);	
		}
		else
		{
			my @line = split(/\t/,$_);
			my @var = split(/\;/,$line[0]);
			my $normal_freq = $line[2];

			for (my $i=3;$i<@line;$i++)
			{
				my @metrics_values = split(/\:/,$line[$i]);
				
				#next if($line[$i]==0);
				next if ($metrics_values[0] == 0);		

				my $out_file =  './VCF/multiSNV/' . $case . '/' . $names[$i-2] . '_multiSNV_sample.vcf';
				if(!(-e $out_file))
				{
					open(O,">>$out_file");
					print O "\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$names[$i-2]\t$names[0]\n";
					print O "$var[0]\t$var[1]\t\.\t$var[2]\t$var[3]\t\.\t\.\tIdentified\=multiSNV\tFREQ\:AD\:DP\t$line[$i]\t$line[2]\n";
					close(O);
				}
				else
				{
					open(O,">>$out_file");
					print O "$var[0]\t$var[1]\t\.\t$var[2]\t$var[3]\t\.\t\.\tIdentified\=multiSNV\tFREQ\:AD\:DP\t$line[$i]\t$line[2]\n";
					close(O);
				}
			}
		}
	}
	close(MERGED_MAT);
}
close(D);

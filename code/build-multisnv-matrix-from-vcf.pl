#!/usr/bin/perl

use strict;

sub get_names
{
	my $string_names = shift;

        my @file_names = split(/\s+/,$string_names);
        
	my @sample_names;

	for(my $i=0;$i<@file_names;$i++)
        {
		my @name_info = split(/\//,$file_names[$i]);
                $name_info[$#name_info]=~/(\S+)\.filtered\.bam/;
		push @sample_names,$1;
	}
	return \@sample_names;
}

sub print_names
{
	my $string_header = shift;
	my $names = shift;
	my $out_mat = shift;

        my @line = split(/\t/,$string_header);
        
	for(my $i=0;$i<=8;$i++)
        {
        	if($i==0)
                {
                	print $out_mat $line[$i];
                }
                else
                {
                	print $out_mat "\t",$line[$i];
                }
        }
	for(my $i=0;$i<@{$names};$i++)
	{
		print $out_mat "\t",@{$names}[$i];
	}
	print $out_mat "\n";
}

my %bp2index=(A=>'0',C=>'1',G=>'2',T=>'3');

open(D,"../list.patients") || die "$!\n";

while(<D>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	
	my $case = $line[0];
	print STDERR $case,"\n";

	my @matrices = `ls ./VCF/multiSNV/$case/\*_multiSNV.vcf`;

	my $printed_header=0;
	
	my $out_matrix_file = "\.\/VCF\/multiSNV\/$case\/$case\_multiSNV_matrix\.txt";

	if(-e $out_matrix_file)
	{
		print STDERR "Removing old version $out_matrix_file\n";
		system("rm $out_matrix_file");
	}

	open my $out_matrix_fh, '>>', $out_matrix_file;


	for(my $k=0;$k<@matrices;$k++)
	{

		chomp($matrices[$k]);
		open(VCF,$matrices[$k]) || die "$!\n";

		my $out_file = $matrices[$k];

		$out_file=~s/\.vcf/\_filtered_renamed.vcf/;

		open my $out_mat_fh, '>', $out_file;


		my @sample_names;

		while(<VCF>)
		{
			chomp($_);
			
			if ($_=~/^\#\#Input\ bam\ files\ \=\ (.+)$/)
			{
				my $string_names = $1;
				my $names = get_names($string_names);
				@sample_names = @{$names};

			}
			if($_=~/^\#/ & $_!~/^\#CHROM/)
			{
				print $out_mat_fh $_,"\n";
			}
			elsif($_=~/^\#/ & $_=~/^\#CHROM/)
			{
				print_names($_,\@sample_names,$out_mat_fh);

			}
			else
			{
				my @line = split(/\t/,$_);

				next if (($line[6] ne 'PASS') & ($line[6] ne 'NORMAL_CONTAMINATION'));

				print $out_mat_fh $_,"\n";

				my @alts = split(/\,/,$line[4]);
				my @metrics_name = split(/\:/,$line[8]);

				#get AF
				my $index_bcount = -1;
				my $index_dp = -1;
				for(my $i=0;$i<@metrics_name;$i++)
				{
					if($metrics_name[$i] eq 'BCOUNT')
					{
						$index_bcount = $i;
					}
					if($metrics_name[$i] eq 'DP')
					{
						$index_dp = $i;
					}               
				}
			
				if($index_bcount == -1 | $index_dp == -1)
				{
					print STDERR "WARN: index_bcount or index_bp not found\n";
					last;
				}


				my $coord;
		
				if($printed_header == 0)
				{		
					print $out_matrix_fh "var\tflag";
					for(my $i=0;$i<@sample_names;$i++)
					{
						print $out_matrix_fh "\t$sample_names[$i]";
					}
					print $out_matrix_fh "\n";
					$printed_header=1;
				}

				#handles coordinates with more than one variant called 
				for(my $i=0;$i<@alts;$i++)
				{
					$coord = "$line[0]\;$line[1]\;$line[3]\;$alts[$i]";
					my $to_print = "";
					for(my $j=9;$j<@line;$j++)
					{
						my @metrics = split(/\:/,$line[$j]);
						my @alt_count = split(/\,/,$metrics[$index_bcount]);
						#print STDERR "$line[0]\;$line[1]\;$line[3]\;$alts[$i]\t$alt_count[$bp2index{$alts[$i]}]\t$metrics[$index_dp]\n";
						my $alt_freq = $alt_count[$bp2index{$alts[$i]}]/$metrics[$index_dp];
						$to_print.= "\t$alt_freq\:$alt_count[$bp2index{$alts[$i]}]\:$metrics[$index_dp]";
					}
					print $out_matrix_fh "$coord\t$line[6]$to_print\n";
				}
			}
		}
		close(VCF);
	}
}
close(D);

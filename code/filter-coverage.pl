#!/usr/bin/perl

use strict;

open(MAT,$ARGV[0]) || die "$!\n";
open(COV,$ARGV[1]) || die "$!\n";


my $th = $ARGV[2];

my %keep=();

while(<COV>)
{
	chomp($_);

	next if ($_!~/\t\d+\t/);

	my @line = split(/\t/,$_);

	my $to_skip = 0;
	
	for(my $i=2;$i<@line;$i++)
	{
		next if ($line[$i]>=$th);
		$to_skip = 1;
	}

	if ($to_skip == 0)
	{
		$keep{"$line[0]\t$line[1]"}++;
	}
}
close(COV);

while(<MAT>)
{
	chomp($_);
	if($_=~/sample/) 
	{
		print $_,"\n";
		next;
	}
	my @line = split(/\t/,$_);
	next if (!defined $keep{"$line[2]\t$line[3]"});
	print $_,"\n";
}
close(MAT);

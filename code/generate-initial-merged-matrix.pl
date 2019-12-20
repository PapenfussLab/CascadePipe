#!/usr/bin/perl

use strict;

my $case = $ARGV[0];
my $germline_name  = $ARGV[1];
my $gender = $ARGV[2]; 
my $cov_th = $ARGV[3];
my $platform = $ARGV[4];

my @files = `ls ./VCF/pmcc_pipeline/$case/\*comb.vcf`;


sub extract_variants_multisnv
{
	my $multisnv_file = shift;
	my $germline_name = shift;
	my $tumor_name = shift;
	my $gender = shift;
	my $cov_th = shift;

	my %multisnv_vars;
 	open(P,$multisnv_file) || die "MultiSNV: $!\n";

	my $index_tumor;
	my $index_normal;

        while(<P>)
        {
                chomp($_);
		if($_=~/^\#/)
		{
                	if ($_=~/^\#CHROM/)
			{
				($index_normal,$index_tumor) = check_index_tumor_vcf($_,$germline_name,$tumor_name);
			}
		}
		else
		{
                	my @line = split(/\t/,$_);
			next if ($line[3]=~/N/ | $line[4]=~/N/);
			my @values = split(/\:/,$line[$index_tumor]);
			my @values_normal = split(/\:/,$line[$index_normal]);

			next if ($values_normal[0]>0.05);
			next if ($values_normal[2]<$cov_th);
			next if ($values[2]<$cov_th);
			next if (($line[0] eq 'Y') & ($gender eq 'F'));

                	$multisnv_vars{"$line[0]\;$line[1]\;$line[3]\;$line[4]"}="multisnv\t$values[0]\t$values[2]\t$values[1]\t$values_normal[0]\t$values_normal[2]\t$values_normal[1]";
        	}
	}
        close(P);
	print STDERR "MultiSNV: ",$index_normal,"\t",$index_tumor,"\n";
	return \%multisnv_vars;
	
}

sub check_index_tumor_vcf
{
	my $header = shift;
	my $germline_name=shift;
	my $tumor_name=shift;

	my @header_split = split(/\t/,$header);

 	my $index_tumor = -1;	
	my $index_germline = -1;

	for(my $i=0;$i<@header_split;$i++)
	{
		if($header_split[$i] eq $germline_name)
		{
			$index_germline = $i;
		}
		if($header_split[$i] eq $tumor_name)
		{
			$index_tumor = $i;
		}
	}

	if($index_tumor == -1 | $index_germline == -1)
	{
		print STDERR "Warn: Normal or Tumour sample not found in VCF\n";
		last;
	}

	return ($index_germline,$index_tumor);

}


sub extract_variants_mpileup_snp
{
	my $mpileup_snp_file = shift;
	my $germline_name = shift;
	my $tumor_name = shift;
	my $gender = shift;
	my $cov_th = shift;

	my %mpileup_snp_vars;
	
	open(VS_SNP,$mpileup_snp_file) || die "$!\n";

	my ($index_tumor,$index_normal);

        while(<VS_SNP>)
        {
                chomp($_);
                my @line = split(/\t/,$_);

		#need to distinguish relaxed varscan calls generated for WES and WGS
		if($_=~/^#CHROM/ & scalar(@line)==11)
		{
			($index_normal,$index_tumor) = (9,10);
		}
		elsif($_=~/^#CHROM/ & scalar(@line)>11)
		{
			($index_normal,$index_tumor) = check_index_tumor_vcf($_,$germline_name,$tumor_name);
		}

                next if ($_=~/^#/);

		#basic filter and rather unreliable, according to Dan Koboldt
		#https://sourceforge.net/p/varscan/discussion/1073559/thread/961ed493/
		#next if ($line[6]!~/PASS/);		#variant did not pass (strand bias or indel error)
		next if($line[$index_tumor] =~ /^\.\/\./);        #tumor doesnt have info

#COMMENTED 07102019, AS IT MAY BE CAUSING FALSE NEGATIVES AT LOW FRACTIONS
#		next if($line[$index_tumor] =~ /^0\/0/);          #tumor has reference homozygous
		next if ($line[4]=~/\,|\+|\-/);		#two variant alleles but varscan output reports one set of values, or indels (handled separately) - skip	
		next if ($line[3]=~/N/ | $line[4]=~/N/);
		next if (($line[0] eq 'Y') & ($gender eq 'F'));

		my @metrics = split(/\:/,$line[8]);
		my $index_freq  = -1;
		my $index_dp  = -1;
		my $index_ad = -1;

		for(my $j=0;$j<@metrics;$j++)
		{
			if ($metrics[$j] eq "FREQ")
			{
				$index_freq  = $j;
			}
			if ($metrics[$j] eq "DP")
			{
				$index_dp  = $j;
			}
			if ($metrics[$j] eq "AD")
			{
				$index_ad  = $j;
			}
		}
		
		last if (($index_freq == -1) | ($index_dp == -1)  | ($index_ad == -1));


		#next if($line[10] !~ /^(0|1)\/(0|1)/);  #tumor doesnt have mutation
                my @values=split(/\:/,$line[$index_tumor]);
                my $allelic_freq = $values[$index_freq];
                $allelic_freq=~s/\%//;

                $allelic_freq = $allelic_freq/100;

		my @values_normal = split(/\:/,$line[$index_normal]);
		my $normal_freq = $values_normal[$index_freq];
		$normal_freq=~s/\%//;
		$normal_freq= $normal_freq/100;

		next if ($normal_freq>0.05);
		next if ($values[$index_dp]<$cov_th);
		next if ($values_normal[$index_dp]<$cov_th);

                $mpileup_snp_vars{"$line[0]\;$line[1]\;$line[3]\;$line[4]"} = "VS_SNV\t$allelic_freq\t$values[$index_dp]\t$values[$index_ad]\t$normal_freq\t$values_normal[$index_dp]\t$values_normal[$index_ad]";
        }
	print STDERR "VS_SNP: ",$index_normal,"\t",$index_tumor,"\n";
        close(VS_SNP);
	return \%mpileup_snp_vars;

}

sub extract_variants_mpileup_indel
{
	my $mpileup_indel_file = shift;
	my $germline_name = shift;
	my $tumor_name = shift;
	my $gender = shift;
	my $cov_th = shift;

	my %mpileup_indel_vars;
	
	open(VS_INDEL,$mpileup_indel_file) || die "$!\n";

	my ($index_normal,$index_tumor);

        while(<VS_INDEL>)
        {
                chomp($_);
                my @line = split(/\t/,$_);

		#need to distinguish relaxed varscan calls generated for WES and WGS
		if($_=~/^#CHROM/ & scalar(@line)==11)
		{
			($index_normal,$index_tumor) = (9,10);
		}
		elsif($_=~/^#CHROM/ & scalar(@line)>11)
		{
			($index_normal,$index_tumor) = check_index_tumor_vcf($_,$germline_name,$tumor_name);
		}


                next if ($_=~/^#/);


		#next if ($line[6]!~/PASS/);		#variant did not pass (strand bias or indel error)
#COMMENTED 07102019, AS IT MAY BE CAUSING FALSE NEGATIVES AT LOW FRACTIONS
#		next if($line[$index_tumor] =~ /^0\/0/);          #tumor has reference homozygous
		next if($line[$index_tumor] =~ /^\.\/\./);        #tumor doesnt have info
		next if ($line[4]=~/\,/);               #two variant alleles but varscan output reports one set of values - skip 
		next if ($line[3]=~/N/ | $line[4]=~/N/);
		next if (($line[0] eq 'Y') & ($gender eq 'F'));
		
		my @metrics = split(/\:/,$line[8]);
		my $index_freq  = -1;
		my $index_dp  = -1;
		my $index_ad = -1;

		for(my $j=0;$j<@metrics;$j++)
		{
			if ($metrics[$j] eq "FREQ")
			{
				$index_freq  = $j;
			}
			if ($metrics[$j] eq "DP")
			{
				$index_dp  = $j;
			}
			if ($metrics[$j] eq "AD")
			{
				$index_ad  = $j;
			}
		}
		
		last if (($index_freq == -1) | ($index_dp == -1)  | ($index_ad == -1));

                my @values=split(/\:/,$line[$index_tumor]);
                my $allelic_freq = $values[$index_freq];
                $allelic_freq=~s/\%//g;

                if($line[4]=~/^\+/)
                {
			$line[4]=~s/^\+/$line[3]/;
		}
		elsif($line[4]=~/^\-(\S+)/)
		{
			my $del = $1;
			$line[4] = $line[3];
			$line[3] = $line[3] . $del;
		}
		else
		{
			next;
		}

		$allelic_freq = $allelic_freq/100;
		
		my @values_normal = split(/\:/,$line[$index_normal]);
		my $normal_freq = $values_normal[$index_freq];
		$normal_freq=~s/\%//;
		$normal_freq= $normal_freq/100;

		next if ($normal_freq>0.05);	
		next if ($values[$index_dp]<$cov_th);
		next if ($values_normal[$index_dp]<$cov_th);
	
		$mpileup_indel_vars{"$line[0]\;$line[1]\;$line[3]\;$line[4]"} = "VS_INDEL\t$allelic_freq\t$values[$index_dp]\t$values[$index_ad]\t$normal_freq\t$values_normal[$index_dp]\t$values_normal[$index_ad]";
        }
	print STDERR "VS_INDEL: ",$index_normal,"\t",$index_tumor,"\n";
        close(VS_INDEL);
	return \%mpileup_indel_vars;

}

sub extract_variants_pmcc
{
        my $mpileup_pmcc_file = shift;
	my $germline_name = shift;
	my $tumor_name = shift;
	my $gender = shift;
	my $cov_th = shift;

        my %mpileup_pmcc_vars;

        open(PMCC,$mpileup_pmcc_file) || die "$!\n";

	my %variants_pmcc=();

        my $index_tumor;
        my $index_normal;

        while(<PMCC>)
        {
                chomp($_);
                if ($_=~/^#CHROM/)
                {
			($index_normal,$index_tumor) = check_index_tumor_vcf($_,$germline_name,$tumor_name);
		}
                next if ($_=~/^#/);

                my @line = split(/\t/,$_);

                my $chrom = $line[0];
                my $coord = $line[1];

		next if($line[$index_tumor] =~ /^\.\/\./);
                next if($line[$index_tumor] =~ /^0\/0/);
		next if (($line[0] eq 'Y') & ($gender eq 'F'));
		next if ($line[3]=~/N/ | $line[4]=~/N/);
	        
		my @metrics = split(/\:/,$line[8]);
		my $index_pmcfreq  = -1;
		my $index_pmcdp  = -1;
		my $index_pmcad = -1;
		my $index_pmcrd = -1;

		for(my $j=0;$j<@metrics;$j++)
		{
			if ($metrics[$j] eq "PMCFREQ")
			{
				$index_pmcfreq  = $j;
			}
			if ($metrics[$j] eq "PMCDP")
			{
				$index_pmcdp  = $j;
			}
			if ($metrics[$j] eq "PMCAD")
			{
				$index_pmcad  = $j;
			}
			if ($metrics[$j] eq "PMCRD")
			{
				$index_pmcrd  = $j;
			}
		}
		
		last if (($index_pmcfreq == -1) | ($index_pmcdp == -1)  | ($index_pmcad == -1) | ($index_pmcrd == -1));

                my @values=split(/\:/,$line[$index_tumor]);
                my $allelic_freq = $values[$index_pmcfreq];
                my $ref_depth = $values[$index_pmcrd];
                my $alt_depth = $values[$index_pmcad];
                my $total_depth = $values[$index_pmcdp];

		my @values_normal = split(/\:/,$line[$index_normal]);

		my $pmc_freq_normal = $values_normal[$index_pmcfreq];

		$line[7]=~s/Identified\=//;
		
		next if ($total_depth < $cov_th);
		next if ($values_normal[$index_pmcdp] < $cov_th);

                my (@multi_alleles,@multi_freqs_tum,@multi_alt_depth_tum,@multi_freqs_normal,@multi_alt_depth_normal);
                if($line[4] =~/\,/)
                {
                        @multi_freqs_tum=split(/\|/,$allelic_freq);
                        @multi_alleles=split(/\,/,$line[4]);
                        @multi_alt_depth_tum = split(/\|/,$alt_depth);
			@multi_freqs_normal = split(/\|/,$pmc_freq_normal);
			@multi_alt_depth_normal = split(/\|/,$values_normal[$index_pmcad]);
			

                        for(my $k=0;$k<@multi_alleles;$k++)
                        {
				next if ($multi_freqs_normal[$k]>0.05);

				next if (($ref_depth + $multi_alt_depth_tum[$k]) < (2*$total_depth/3));
                                next if ((length($line[3])!=1 | length($multi_alleles[$k])!=1) & $multi_freqs_tum[$k]<0.2);
                                $variants_pmcc{"$line[0]\;$line[1]\;$line[3]\;$multi_alleles[$k]"} = "$line[7]\t$multi_freqs_tum[$k]\t$total_depth\t$multi_alt_depth_tum[$k]\t$multi_freqs_normal[$k]\t$values_normal[$index_pmcdp]\t$multi_alt_depth_normal[$k]";
			}
		}
                else
                {
			next if ($pmc_freq_normal>0.05);
                	#this to deal with SNVs overlapping with InDels, which has a high total_depth, but a very low ref_depth and alt_depth and likely FP Indels
                	next if (($ref_depth + $alt_depth) < (2*$total_depth/3));
                        next if ((length($line[3])!=1 | length($line[4])!=1) & $allelic_freq<0.2);
                        $variants_pmcc{"$line[0]\;$line[1]\;$line[3]\;$line[4]"} = "$line[7]\t$allelic_freq\t$total_depth\t$alt_depth\t$pmc_freq_normal\t$values_normal[$index_pmcdp]\t$values_normal[$index_pmcad]";
		}
	}
	print STDERR "PMCC: ",$index_normal,"\t",$index_tumor,"\n";
	close(PMCC);
	return \%variants_pmcc;
}
	

print "sample\tvariant\tchromosome\tcoordinate\tref\talt\tcategory\tvariant_caller\ttumour_af\ttumour_dp\ttumour_ad\tgermline_af\tgermline_dp\tgermline_ad\n";

for(my $i=0;$i<@files;$i++)
{
	chomp($files[$i]);
	my @name_info = split(/\//,$files[$i]);

	$name_info[$#name_info]=~/(\S+)_merged/;
	my $sample = $1;
	
	my $multisnv = './VCF/multiSNV/' . $case . '/' . $sample . '_multiSNV.vcf';

	my ($mpileup_snp,$mpileup_indel);

	if($platform eq 'WES')
	{
		$mpileup_snp = './VCF/mpileup/' . $case . '/' . $sample . '.mpileup.snp.vcf';
		$mpileup_indel = './VCF/mpileup/' . $case . '/' . $sample . '.mpileup.indel.vcf';
	}
	elsif($platform eq 'WGS')
	{
		$mpileup_snp = './VCF/mpileup/' . $case . '/' . $case . '.mpileup.snp.vcf';
		$mpileup_indel = './VCF/mpileup/' . $case . '/' . $case . '.mpileup.indel.vcf';
	}
	else
	{
		print STDERR "Platform needs to be WES or WGS (due to way relaxed varscan was called)\n";
		last;
	}

	print STDERR $multisnv,"\n";

	my $pmcc_vars = extract_variants_pmcc($files[$i],$germline_name,$sample,$gender,$cov_th);
	my $multisnv_vars = extract_variants_multisnv($multisnv,$germline_name,$sample,$gender,$cov_th);
	my $mpileup_snp_vars = extract_variants_mpileup_snp($mpileup_snp,$germline_name,$sample,$gender,$cov_th);
	my $mpileup_indel_vars = extract_variants_mpileup_indel($mpileup_indel,$germline_name,$sample,$gender,$cov_th);


	for my $key (keys %$pmcc_vars)
	{
		my $pmcc_info = $pmcc_vars->{$key};
		my @var_info = split(/\;/,$key);
	
		#ADDED 07102019: WILL INCLUDE MUTECT CALLS AS CATEGORY 2 SINCE IT IS SENSITIVE FOR LOW PURITY SAMPLES
		if($pmcc_info=~/^mutect\t/ & (!defined $multisnv_vars->{$key}))
		{
			 print $sample,"\t",$key,"\t$var_info[0]\t$var_info[1]\t$var_info[2]\t$var_info[3]\t2\t",$pmcc_info,"\n";
			
		}
	
		#variant is supported by only one variant caller
		next if ($pmcc_info!~/\-/ & (!defined $multisnv_vars->{$key}));


		if(defined $multisnv_vars->{$key})
		{
			print $sample,"\t",$key,"\t$var_info[0]\t$var_info[1]\t$var_info[2]\t$var_info[3]\t1\tmultisnv\-",$pmcc_info,"\n";
		}
		else
		{
			 print $sample,"\t",$key,"\t$var_info[0]\t$var_info[1]\t$var_info[2]\t$var_info[3]\t1\t",$pmcc_info,"\n";
		}
	}

	for my $key (keys %$multisnv_vars)
	{
		my $multisnv_info = $multisnv_vars->{$key};
		next if (defined $pmcc_vars->{$key});
		next if (!defined $mpileup_snp_vars->{$key});

		my @var_info = split(/\;/,$key);

		print $sample,"\t",$key,"\t$var_info[0]\t$var_info[1]\t$var_info[2]\t$var_info[3]\t2\tVS_SNV\-",$multisnv_info,"\n";
	}

	for my $key (keys %$mpileup_indel_vars)
	{
		next if (defined $pmcc_vars->{$key});
		my $mpileup_indel_info = $mpileup_indel_vars->{$key};
		my @var_info = split(/\;/,$key);
		print $sample,"\t",$key,"\t$var_info[0]\t$var_info[1]\t$var_info[2]\t$var_info[3]\t2\t",$mpileup_indel_info,"\n";
	}

	#from here, need to work with these four hash to get as output:
	#(1) sample, supporting variant callers, allele fraction tumour, allele fraction normal, variant read count, total read count
	#(2) then a separate script will build the matrix from the output

	#say $_, " => ", $result->{$_} for keys %$result;


}


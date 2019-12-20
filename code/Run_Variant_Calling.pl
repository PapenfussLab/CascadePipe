#!/usr/bin/perl

#Created by Jason Ellul, April 2013


use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use Pipeline::Utils;
use File::Slurp;
use YAML::XS;
use IPC::System::Simple qw(capture);
use Pipeline::VariantCalling::Utils;
use Text::ASCIITable;
$| = 1;

# Grab and set all options
my %OPTIONS = ( f => 0, v => "ug" );
getopts('ab:cde:f:ghin:o:p:q:rsv:w:x', \%OPTIONS);

die qq(
Usage:   Run_Variant_Calling.pl [OPTIONS] <sample sheet> <genome/amplicon set> <analysis type> <email>

OPTIONS:		-a		Display Analysis Types Available
			-b	STR	Bed file
			-c		Output GVCF file if using GATK HaplotypeCaller
			-d		Mark Duplicates
			-e	STR	Ensembl VEP plugins and their config options to run, plugins must be seperated by a colon e.g. Plugin1,/path/to/config/Plugin1/config,b:Plugin2
					Currently only Genome_Trax is installed for all users as an additional plugin. See http://asia.ensembl.org/info/docs/variation/vep/vep_script.html for how
					to install a plugin in your home directory
			-f	STR	Set the familial output folder
			-g		Keep germline variants when performing Somatic analysis
			-h		Remove LOH calls from somatic variants
			-i		Perform indel realignment
			-n	STR	Name to prefix output (only used for familial)
			-o	STR	Overide default genome build (WARNING this will use an external database which may be extremely slow)
			-p	STR	Directory to place the pbs files [<path to sample sheet>]
			-q	STR	Queue to submit jobs to
			-r		Perfrom GATK Base Quality Recalibration
			-s		Sort the variant file by Transcript and Consequence
			-v	STR	List of variant callers to use for the analysis (available: hap - GATK HaplotypeCaller, ss - Somatic Sniper, jsnvm jointSNVMix, mutect - MuTect, id - IndelGenotyper, pindel - Pindel, vs - Varscan, ug - GATK Unified Genotyper) [$OPTIONS{v}]
			-w	INT	Walltime scale factor
			-x		Create pbs script but don't run it

Sample Sheet format:
Familial:	bam_file
Somatic:	output_base_lcation	normal_bam_file	tumor_bam_file
Other:		output_base_lcation	bam_file

) unless @ARGV == 4 || defined $OPTIONS{a};

# version info
my $version = 0.3;
my $Script = 'Run_Variant_Calling';

my $default_db = "78";
my $db_version = defined $OPTIONS{o} ? $OPTIONS{o} : $default_db;

# files
my $References = Load scalar read_file("/data/reference/indexes/Index_map.txt");
my $BedFiles = Load scalar read_file("/data/reference/bed_files/Bedfile_map.txt");
my $G1000_omni = "/data/databases/Gatk_Reseource_Bundle/1000G_omni2.5.b37.vcf";
my $G1000_indels = "/data/databases/Gatk_Reseource_Bundle/1000G_phase1.indels.b37.vcf";
my $dbsnp = "/data/databases/Gatk_Reseource_Bundle/dbsnp_137.b37.vcf";
my $mouse_dbsnp = "/data/databases/dbSNP/mouse_dbsnp137/Mus_musculus.vcf";
my $fly_dbsnp = "/data/databases/dbSNP/fruitfly_dpgp_1/Drosophila_melanogaster.vcf";
my $hapmap = "/data/databases/Gatk_Reseource_Bundle/hapmap_3.3.b37.vcf";
my $mills = "/data/databases/Gatk_Reseource_Bundle/Mills_and_1000G_gold_standard.indels.b37.vcf";
my $cosmic = "/data/databases/Cosmic/Cosmic_v54/hg19_cosmic_v54_120711.vcf";
my %variant_caller_names = (hap => "GATK HaplotypeCaller", ss => "Somatic Sniper", jsnvm => "jointSNVMix", mutect => "MuTect", id => "IndelGenotyper", pindel => "Pindel", vs => "Varscan", ug => "GATK Unified Genotyper");
my %germline_callers = ( vs => 1, pindel => 1, ug => 1, hap => 1 );
my %Analysis_Type;
foreach my $key (keys %{ $BedFiles }) {
	my $new = lc($key);
	$new =~ s/ /_/g;
	$Analysis_Type{$new} = $key;
}

if(defined $OPTIONS{o}) {
	if($OPTIONS{o} == 67 || $OPTIONS{o} == 69 || $OPTIONS{o} == 73) {
		warn "WARNING overiding the genome build will use an external database to annotate your data which may be extremely slow. Caches exists for builds 67, 69 and 73.";
	} else {
		die "Genome build not known, caches exist only for builds 67, 69 and 73.";
	}
}

# walltime scaling
my $walltime = defined $OPTIONS{w}? $OPTIONS{w} : 1;

# queue
my $queue = defined $OPTIONS{q}? " -q $OPTIONS{q}" : "";

# display analysis types
if(defined $OPTIONS{a}) {
	foreach my $new (sort keys %Analysis_Type) {
		my $anlysis = $Analysis_Type{$new};
		if(defined $$References{$$BedFiles{$anlysis}{aligner}}) {
			print "$new:\n";
			if(defined $$BedFiles{$anlysis}{var}) {
				print "\t$_\n" foreach sort keys %{ $$BedFiles{$anlysis}{var} };
			} else {
				print "\t$_\n" foreach sort keys %{ $$References{$$BedFiles{$anlysis}{aligner}} };
			}
		}
	}
	exit;
}

# grab inputs
my $samplesheet = shift;
my $genome = shift;
my $type = shift;
my $email = shift;
$type = $Analysis_Type{$type} if defined $Analysis_Type{$type};

# check type
die("Analysis type: $type not known, exiting!\n") unless defined $$BedFiles{$type}{aligner};
my $aligner = $$BedFiles{$type}{aligner};

# check Genome
$genome = $$References{expand}{$genome} if defined $$References{expand}{$genome} && $aligner ne "amplicon" && $aligner ne "amplicon_halo";
die("Genome: $genome does not exist\n") unless defined $$References{$aligner}{$genome};
my $ref = $$References{fasta}{$genome};

# grab bed file
$OPTIONS{b} = $$BedFiles{$type}{var}{$genome} unless defined $OPTIONS{b} || !defined $$BedFiles{$type}{var}{$genome};
$OPTIONS{b} = "" unless defined $OPTIONS{b};

# for amplicon now get the actual genome rather than the amplicon set
$genome = $$References{expand}{$genome} if defined $$References{expand}{$genome};

# determine species
my $species = $$References{species}{$genome};

# print what we are using for variant calling
print "Using for variant calling:\n\tGenome:      $genome\n\tReference:   $ref\n\tSpecies:     $species\n\tBedfile:     $OPTIONS{b}\n";

# options
my %tumor;
my @tum_ord;
my @norm_ord;
$OPTIONS{c} = defined $OPTIONS{c}? 1:0;
$OPTIONS{d} = defined $OPTIONS{d}? 1:0;
$OPTIONS{g} = defined $OPTIONS{g}? 1:0;
$OPTIONS{h} = defined $OPTIONS{h}? 1:0;
$OPTIONS{i} = defined $OPTIONS{i}? 1:0;
$OPTIONS{r} = defined $OPTIONS{r}? 1:0;
$OPTIONS{s} = defined $OPTIONS{s}? 1:0;

my $GATKversion;
if($OPTIONS{i} || $OPTIONS{r}) {
	$GATKversion = `GenomeAnalysisTK -h | grep \"Genome Analysis Toolkit (GATK)\"`;
	chomp $GATKversion;
	$GATKversion =~ s/The Genome Analysis Toolkit \(GATK\) (.*), Compiled .*/$1/;
}

my $PicardVerion;
if($OPTIONS{d}) {
	$PicardVerion = `BuildBamIndex.sh --version 2>&1`;
	chomp $PicardVerion;
}

my $err;
# check parameters for know data types
if($$BedFiles{$type}{aligner} eq "amplicon" || $$BedFiles{$type}{aligner} eq "amplicon_halo") {
	if($OPTIONS{d}) {
		$err = "Warning: Mark Duplicates should not be used for amplicon data, over riding and setting Mark Duplicates to off!\n";
		$OPTIONS{d} = 0;
	}
	if($OPTIONS{r}) {
		$err .= "Warning: GATK Base Quality Recalibration should not be used for amplicon data, over riding and setting GATK Base Quality Recalibration to off!\n";
		$OPTIONS{r} = 0;
	}
} elsif($$BedFiles{$type}{aligner} eq "haloplex") {
	if($OPTIONS{d}) {
		$err = "Warning: Mark Duplicates should not be used for haloplex data, over riding and setting Mark Duplicates to off!\n";
		$OPTIONS{d} = 0;
	}
	if($OPTIONS{r}) {
		print "Warning: GATK Base Quality Recalibration should not be used for haloplex data, are you sure you wish to proceed?\n";
		my $quit = 0;

		until ($quit) {
    		print "Enter Y|N|Q: ";
    		chomp(my $input = <STDIN>);
    
    		if ($input =~ /^[Yy]?$/i) {      # Match Yy or blank
       			$err .= "Warning: Applying GATK Base Quality Recalibration to haloplex data\n";
       			$quit = 1;
    		} elsif ($input =~ /^[Nn]$/i) {  # Match Nn
        		print "Continuing without GATK Base Quality Recalibration\n";
        		$OPTIONS{r} = 0;
        		$quit = 1;
    		} elsif ($input =~ /^[Qq]$/i) {  # Match Qq
        		print "Quitting.\n";
        		exit;
    		} else {
        		print "Enter Y|N|Q: ";
    		}
		}		
	}
}

# read the sample sheet file
open(IN, $samplesheet) or die("$samplesheet: $!\n");
my $analysis_type = "";
while(<IN>) {
	chomp;

	# skip any empty lines
	next if $_ eq "";

	# split line
	my @samples = split("\t", $_);
	$analysis_type = "familial" if $analysis_type eq "" && @samples == 1;
	$analysis_type = "standard" if $analysis_type eq "" && @samples == 2;
	$analysis_type = "tumor_norm" if $analysis_type eq "" && @samples == 3;

	# check output dir is specified for familial analysis
	die("For Familial analyses the output dir must be specified using the -f parameter!\n") unless $OPTIONS{f} || $analysis_type ne "familial";

	# check all samples are of the same type
	die("Error: All samples need to be the same, for familial only the bam needs to be listed, output base lcation must be specified using the -f parameter!\n") if @samples != 1 && $analysis_type eq "familial";
	die("Error: All samples need to be the same, for standard variant calling the file should contain 2 columns (output base lcation and then sample bam file)!\n") if @samples != 2 && $analysis_type eq "standard";
	die("Error: All samples need to be the same, for tumor normal variant calling the file should contain 3 columns (output base lcation, normal bam file and then tumor bam file)!\n") if @samples != 3 && $analysis_type eq "tumor_norm";

	# check bamfiles
	die("Error: Sample $samples[0] do not exist, please check sample file\n") unless $analysis_type ne "familial" || -f $samples[0];
	die("Error: Sample $samples[1] do not exist, please check sample file\n") unless $analysis_type eq "familial" || -f $samples[1];
	die("Error: Sample $samples[2] do not exist, please check sample file\n") unless $analysis_type ne "tumor_norm" || -f $samples[2];

	# get full path name of directories and files
	my $index;
	for($index = 0; $index <= $#samples; $index++) {
		$samples[$index] = capture("readlink -f $samples[$index]");
		chomp $samples[$index];
	}
	$index--;

	# extract sample name from bam file and save folder and check if sample info already loaded
	my $name = bam_sample_name($samples[$index]);
	die("Error: The sample $name is specified more than once, $samples[$index] && $tumor{$name}{bam_file}!\n") if defined $tumor{$name}{bam_file};
	$tumor{$name}{bam_file} = $samples[$index];

	# we have tumor and normal samples
	if($analysis_type eq "tumor_norm") {
		# extract sample name from bam file and save folder
		my $norm = bam_sample_name($samples[1]);
		$tumor{$name}{Norm} = $norm;
		$tumor{$name}{norm_bam_file} = $samples[1];
	}

	if($analysis_type ne "familial") {
		# set directories for non familial analyses
		$tumor{$name}{Out} = $samples[0];
		$tumor{$name}{Bam} = "$samples[0]/Bam";
		$tumor{$name}{Var} = "$samples[0]/Variants";
		$tumor{$name}{VEP} = "$samples[0]/VEP";
		$tumor{$name}{VCF} = "$samples[0]/VCF";
		$tumor{$name}{Log} = "$samples[0]/Logs";
		$tumor{$name}{TMP} = "$samples[0]/TMP";
	} else {
		# set directories for familial analyses
		mkdir $OPTIONS{f} or die("Unable to create directory $OPTIONS{f}: $!\n") unless -d $OPTIONS{f};
		$OPTIONS{f} = capture("readlink -f $OPTIONS{f}");
		chomp $OPTIONS{f};

		$tumor{$name}{Out} = $OPTIONS{f};
		$tumor{$name}{Bam} = "$OPTIONS{f}/Bam";
		$tumor{$name}{Var} = "$OPTIONS{f}/Variants";
		$tumor{$name}{VEP} = "$OPTIONS{f}/VEP";
		$tumor{$name}{VCF} = "$OPTIONS{f}/VCF";
		$tumor{$name}{Log} = "$OPTIONS{f}/Logs";
		$tumor{$name}{TMP} = "$OPTIONS{f}/TMP";
	}
	mkdir $tumor{$name}{Out} or die("Unable to create directory $tumor{$name}{Out}: $!\n") unless -d $tumor{$name}{Out};
	mkdir $tumor{$name}{Bam} or die("Unable to create directory $tumor{$name}{Bam}: $!\n") unless -d $tumor{$name}{Bam};
	mkdir $tumor{$name}{Var} or die("Unable to create directory $tumor{$name}{Var}: $!\n") unless -d $tumor{$name}{Var};
	mkdir $tumor{$name}{VEP} or die("Unable to create directory $tumor{$name}{VEP}: $!\n") unless -d $tumor{$name}{VEP};
	mkdir $tumor{$name}{VCF} or die("Unable to create directory $tumor{$name}{VCF}: $!\n") unless -d $tumor{$name}{VCF};
	mkdir $tumor{$name}{Log} or die("Unable to create directory $tumor{$name}{Log}: $!\n") unless -d $tumor{$name}{Log};
	mkdir $tumor{$name}{TMP} or die("Unable to create directory $tumor{$name}{TMP}: $!\n") unless -d $tumor{$name}{TMP};
}
close(IN);

# check we have at lease 1 tumor folder, and when somatic that the samples have a matching normal
die("Error: No tumor folders found please respecify\n") unless keys %tumor;

# check varinat callers
my @caller_names;
my @callers = split(",", $OPTIONS{v});
for(my $i = 0; $i <= $#callers; $i++) {
	if(!defined $variant_caller_names{$callers[$i]}) {
		print "Warning: variant caller $callers[$i] is not known!\n";
		splice(@callers, $i, 1);
		$i--;
	} elsif($analysis_type ne "tumor_norm" && !defined $germline_callers{$callers[$i]}) {
		print "Warning: variant caller $callers[$i] is only available for somatic samples!\n";
		splice(@callers, $i, 1);
		$i--;
	} elsif($analysis_type eq "tumor_norm" && $callers[$i] eq 'hap') {
		print "Warning: variant caller $callers[$i] is only available for germline!\n";
		splice(@callers, $i, 1);
		$i--;
	} else {
		push @caller_names, $variant_caller_names{$callers[$i]};
	}
}
die("No valid callers specified!\n") unless @callers;
print "\tCallers:     ".join(", ", @caller_names)."\n\n";
warn("Warning: No somatic variant callers specified even though you have tumor and normal samples, only @caller_names will be used for variant calling!\n") if $analysis_type eq "tumor_norm" && @callers == 1 && $callers[0] eq "ug";

# display the tumor and normal samples found
my $t = Text::ASCIITable->new({ headingText => 'Sample Information' });

my $idx = 1;
if($analysis_type eq "tumor_norm") {
	$t->setCols('Sample', 'Sample Name Tumor', 'Sample Name Normal', 'Out Directory');
	$t->addRow($idx++, $_, $tumor{$_}{Norm}, $tumor{$_}{Out}) foreach sort keys %tumor;
} elsif($analysis_type eq "familial") {
	$t->setCols('Sample', 'Sample Name', 'Bam File');
	$t->addRow($idx++, $_, $tumor{$_}{bam_file}) foreach sort keys %tumor;
	$t->addRowLine();
	$t->addRow('Out Directory', $OPTIONS{f});
} else {
	$t->setCols('Sample', 'Sample Name', 'Out Directory');
	$t->addRow($idx++, $_, $tumor{$_}{Out}) foreach sort keys %tumor;
}
print $t;

# determine the folder for the pbs script
if(!defined $OPTIONS{p}) {
	my ($file, $path) = fileparse($samplesheet);
	$OPTIONS{p} = "$path";
} else {
	mkdir $OPTIONS{p} or die("Cannot create pbs directory $OPTIONS{p}: $!\n") unless -d $OPTIONS{p};
}
$OPTIONS{p} = capture("readlink -f $OPTIONS{p}");
chomp $OPTIONS{p};
print "\nCreating pbs script: $OPTIONS{p}/pbs_submit.sh\n\n";
open(PBS, ">$OPTIONS{p}/pbs_submit.sh") or die("$OPTIONS{p}/pbs_submit.sh: $!\n");

# add the direcory creation to the pbs script just in case we wish to run it again
if($analysis_type eq "familial") {
	print PBS "# Create directories for the output\n";
	print PBS "if [ ! -d \"$OPTIONS{f}/Bam\" ]; then\n\tmkdir \"$OPTIONS{f}/Bam\"\nfi\n";
	print PBS "if [ ! -d \"$OPTIONS{f}/Variants\" ]; then\n\tmkdir \"$OPTIONS{f}/Variants\"\nfi\n";
	print PBS "if [ ! -d \"$OPTIONS{f}/VEP\" ]; then\n\tmkdir \"$OPTIONS{f}/VEP\"\nfi\n";
	print PBS "if [ ! -d \"$OPTIONS{f}/VCF\" ]; then\n\tmkdir \"$OPTIONS{f}/VCF\"\nfi\n";
	print PBS "if [ ! -d \"$OPTIONS{f}/Logs\" ]; then\n\tmkdir \"$OPTIONS{f}/Logs\"\nfi\n";
	print PBS "if [ ! -d \"$OPTIONS{f}/TMP\" ]; then\n\tmkdir \"$OPTIONS{f}/TMP\"\nfi\n\n";
} else {
	foreach my $sample (sort keys %tumor) {
		print PBS "# Create output directories for the sample $sample\n";
		print PBS "if [ ! -d \"$tumor{$sample}{Bam}\" ]; then\n\tmkdir \"$tumor{$sample}{Bam}\"\nfi\n";
		print PBS "if [ ! -d \"$tumor{$sample}{Var}\" ]; then\n\tmkdir \"$tumor{$sample}{Var}\"\nfi\n";
		print PBS "if [ ! -d \"$tumor{$sample}{VEP}\" ]; then\n\tmkdir \"$tumor{$sample}{VEP}\"\nfi\n";
		print PBS "if [ ! -d \"$tumor{$sample}{VCF}\" ]; then\n\tmkdir \"$tumor{$sample}{VCF}\"\nfi\n";
		print PBS "if [ ! -d \"$tumor{$sample}{Log}\" ]; then\n\tmkdir \"$tumor{$sample}{Log}\"\nfi\n";
		print PBS "if [ ! -d \"$tumor{$sample}{TMP}\" ]; then\n\tmkdir \"$tumor{$sample}{TMP}\"\nfi\n\n";
	}
}

if(defined $err) {
	print $err;
}

# for each tumor sample save the original bam location and create symbolic link in tmp folder
symbolic_link(\%tumor, "bam_file");
symbolic_link(\%tumor, "norm_bam_file") if $analysis_type eq "tumor_norm";

# if no bed file set $OPTIONS{b} to undefined
delete $OPTIONS{b} if $OPTIONS{b} eq "";

# mark duplicates
my $move = 0;
my $ext = "";
if($OPTIONS{d}) {
	# mark the duplicate reads
	$ext = "_dups";
	if($analysis_type eq "tumor_norm") {
		markdups(\%tumor, "current_bam_file", "_Tum");
		markdups(\%tumor, "current_norm_bam_file", "_Norm");
	} else {
		markdups(\%tumor, "current_bam_file", "");
	}
}

# merge bam files
my %sep_tumor = %tumor;
if($analysis_type eq "familial") {
	%tumor = ();
	%{ $tumor{Familial} } = %{ $sep_tumor{(keys %sep_tumor)[0]} };
	delete $tumor{Familial}{bam_file};

	# with the familial case we merge all samples into 1 bam file
	print PBS "########### Merging the familial bam files\n";
	my @depend;
	my @bams;
	push @depend, $sep_tumor{$_}{depend} foreach sort keys %sep_tumor;
	push @bams, $sep_tumor{$_}{current_bam_file} foreach sort keys %sep_tumor;
	@{ $tumor{Familial}{merged} } = sort keys %sep_tumor;
	if(defined $OPTIONS{n}) {
		$tumor{Familial}{current_bam_file} = "$tumor{Familial}{TMP}/$OPTIONS{n}_Familial_merged$ext.bam";
	} else {
		$tumor{Familial}{current_bam_file} = "$tumor{Familial}{TMP}/$OPTIONS{n}_Familial_merged$ext.bam";
	}
	$tumor{Familial}{depend} = ":\$MERGE1";
	merge_bams(1, \@depend, \@bams, $tumor{Familial}{current_bam_file}, "$tumor{Familial}{Log}/$OPTIONS{n}_Familial_1_Merge.txt", "$OPTIONS{n}_Familial_1_Merge");

	print PBS "RM1=`qsub_command -m n -N Familial_1_MergeRM -W depend=afterok:\$MERGE1 -j oe -d $OPTIONS{p}$queue 'rm -f @bams'`\n" if $OPTIONS{d};
	print PBS "\n";
} elsif($analysis_type eq "tumor_norm") {
	# with the somatic case merge the tumor and normal pairs
	$idx = 1;
	foreach my $sample (sort keys %tumor) {
		print PBS "########### Merging the tumor normal pair $sample $tumor{$sample}{Norm}\n";
		my @depend = ($tumor{$sample}{depend});
		my @bams = ($tumor{$sample}{current_norm_bam_file}, $tumor{$sample}{current_bam_file});
		@{ $tumor{$sample}{merged} } = ($tumor{$sample}{Norm}, $sample);
		$tumor{$sample}{current_bam_file} = "$tumor{$sample}{TMP}/${sample}_merged$ext.bam";
		$tumor{$sample}{depend} = ":\$MERGE$idx";
		merge_bams($idx, \@depend, \@bams, $tumor{$sample}{current_bam_file}, "$tumor{$sample}{Log}/${sample}_${idx}_Merge.txt", "${sample}_${idx}_Merge");

		print PBS "RM$idx=`qsub_command -m n -N ${sample}_${idx}_MergeRM -W depend=afterok:\$MERGE$idx -j oe -d $OPTIONS{p}$queue 'rm -f @bams'`\n" if $OPTIONS{d};
		$idx++;
		print PBS "\n";
	}
}

# indel realignmnet
if($OPTIONS{i}) {
	$ext .= "_realign";
	indel(\%tumor);
}

# base score recalibration
if($OPTIONS{r}) {
	$ext .= "_recal";
	recal(\%tumor);
}

# run gatk or varscan to detect variants
if(grep {$_ eq "ug"} @callers) {
	ug(\%tumor, $analysis_type);
} else {
	$tumor{$_}{depend_rm} = '' foreach sort keys %tumor;
}

hap(\%tumor) if $analysis_type ne "tumor_norm" && grep {$_ eq "hap"} @callers;
vs_germ(\%tumor) if $analysis_type ne "tumor_norm" && grep {$_ eq "vs"} @callers;
combine_vcf(\%tumor) unless $analysis_type eq "tumor_norm";

# split or move bam file if we have marked dups, realigned or recalibrated
if($analysis_type ne "standard" && $ext ne "") {
	print PBS "########### Splitting the Bam file by sample\n";
	$idx = 1;

	# walltime
	my $wt = int($walltime * 12 * 3600);

	foreach my $sample (sort keys %tumor) {
		print PBS "SPLIT$idx=`qsub_command -m n -N ${sample}_${idx}_SPLIT -l walltime=$wt -W depend=afterok$tumor{$sample}{depend} -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::GATK; SplitSamFile({I => q($tumor{$sample}{current_bam_file}), outputRoot => q($tumor{$sample}{Bam}/Split${ext}_), Log => q($tumor{$sample}{Log}/${sample}_${idx}_SPLIT.txt), R => q($ref)});\\x27'`\n";
		print PBS "RM$idx=`qsub_command -m n -N ${sample}_${idx}_SPLIT_RM -W depend=afterok$tumor{$sample}{depend_rm}:\$SPLIT$idx -j oe -d $OPTIONS{p}$queue 'rm -f $tumor{$sample}{current_bam_file}'`\n";

		$tumor{$sample}{depend} = ":\$SPLIT$idx";
		push @{ $tumor{$sample}{bamstats} },  "$tumor{$sample}{Bam}/Split${ext}_${_}.bam" foreach @{ $tumor{$sample}{merged} };
		$idx++;
	}
	print PBS "\n";
} elsif($ext ne "") {
	# otherwise move the bam files into the sample folders
	print PBS "########### Moving the final bam files into the Bam folder\n";
	$idx = 1;
	foreach my $sample (sort keys %tumor) {
		my $bam = fileparse($tumor{$sample}{current_bam_file});
		my $bai = $tumor{$sample}{current_bam_file};
		$bai =~ s/m$/i/;
		print PBS "MV$idx=`qsub_command -m n -N ${sample}_${idx}_MV -W depend=afterok$tumor{$sample}{depend}$tumor{$sample}{depend_rm} -j oe -d $OPTIONS{p}$queue 'mv $tumor{$sample}{current_bam_file} $tumor{$sample}{Bam}; mv $bai $tumor{$sample}{Bam}'`\n";
		$tumor{$sample}{depend} = ":\$MV$idx";
		push @{ $tumor{$sample}{bamstats} },  "$tumor{$sample}{Bam}/$bam";
		$idx++;
	}
	print PBS "\n";
} elsif($analysis_type eq "tumor_norm") {
	@{ $tumor{$_}{bamstats} } = ($tumor{$_}{norm_bam_file}, $tumor{$_}{bam_file}) foreach sort keys %tumor;
	delete $tumor{$_}{depend} foreach sort keys %tumor;
} elsif($analysis_type eq "familial") {
	push @{ $tumor{Familial}{bamstats} }, $sep_tumor{$_}{bam_file} foreach sort keys %sep_tumor;
	delete $tumor{Familial}{depend};
} else {
	@{ $tumor{$_}{bamstats} } = ($tumor{$_}{bam_file}) foreach sort keys  %tumor;
}

# run pindel if selected
pindel(\%tumor) if grep {$_ eq "pindel"} @callers;

# somatic variant callers
if($analysis_type eq "tumor_norm" && (@callers > 1 || $callers[0] ne "ug")) {
	$idx = 1;
	my $wt = int($walltime * 24 * 3600);
	foreach my $sample (sort keys %tumor) {
		print PBS "########### Running Variant Callers on tumor normal samples $sample\n";
		my $vcf_out = fileparse($tumor{$sample}{current_bam_file});
		$vcf_out =~ s/\.bam$//;
		$vcf_out = "$tumor{$sample}{VCF}/$vcf_out";
		my ($normal, $tumor) = @{ $tumor{$sample}{bamstats} };
		my $germline = "";
		my $depend = "";
		$depend = " -W depend=afterok$tumor{$sample}{depend}" if defined $tumor{$sample}{depend};
		foreach my $caller (@callers) {
			if($caller eq "ss") {
				# somatic sniper
				$germline = ", q(${vcf_out}_germline_ss.vcf)" if $OPTIONS{g};
				my $ss_opt = "somatic_threshold => 40, Log => q($tumor{$sample}{Log}/${sample}_${idx}_SS.txt)";
				$ss_opt .= ", bed => q($OPTIONS{b})" if defined $OPTIONS{b};
				print PBS "SS$idx=`qsub_command -m n -N ${sample}_${idx}_SS -l walltime=$wt$depend -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::SomaticSniper; SomaticSniper(q($normal), q($tumor), [q(${vcf_out}_somatic_loh_ss.vcf)$germline], q($ref), { $ss_opt });\\x27'`\n";

				# remove loh calls into seperate file
				if($OPTIONS{h}) {
					print PBS "SS_loh$idx=`qsub_command -m n -N ${sample}_${idx}_SS_LOH -W depend=afterok:\$SS$idx -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::Utils; filter_tsv_file(q(${vcf_out}_somatic_loh_ss.vcf), q(${vcf_out}_somatic_ss.vcf), { filter => \"[11] !~ \\\\\$.*:3:[^:]*\\\\\$\\\\\$\" });\\x27 > $tumor{$sample}{Log}/${sample}_${idx}_SS_LOH.txt 2>&1'`\n";
					print PBS "SS_Filt$idx=`qsub_command -m n -N ${sample}_${idx}_SS_Filt -W depend=afterok:\$SS_loh$idx -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::Utils; filter_tsv_file(q(${vcf_out}_somatic_loh_ss.vcf), q(${vcf_out}_loh_ss.vcf), { filter => \"[11] =~ \\\\\$.*:3:[^:]*\\\\\$\\\\\$\" });\\x27 > $tumor{$sample}{Log}/${sample}_${idx}_SS_Filt.txt 2>&1'`\n";
					print PBS "RMSS$idx=`qsub_command -m n -N ${sample}_${idx}_SS_RM -W depend=afterok:\$SS_Filt$idx -j oe -d $OPTIONS{p}$queue 'rm -f ${vcf_out}_somatic_loh_ss.vcf'`\n";
					$tumor{$sample}{vcf_comb}{somatic_sniper}{file} = "${vcf_out}_somatic_ss.vcf";
					$tumor{$sample}{vcf_comb}{somatic_sniper}{depend} = ":".insert_version("SS_Filt", "${vcf_out}_somatic_ss.vcf", $sample, "\$SS_Filt$idx", $idx);
				} else {
					$tumor{$sample}{vcf_comb}{somatic_sniper}{file} = "${vcf_out}_somatic_loh_ss.vcf";
					$tumor{$sample}{vcf_comb}{somatic_sniper}{depend} = ":".insert_version("SS", "${vcf_out}_somatic_loh_ss.vcf", $sample, "\$SS$idx", $idx);
				}
			} elsif($caller eq "jsnvm") {
				# jointSNVMix
				$germline = ", \"${vcf_out}_germline_jsvm.vcf\"" if $OPTIONS{g};
				my $jsnvm_opt = "train => 1, Log => q($tumor{$sample}{Log}/${sample}_${idx}_JSNVM.txt)";
				$jsnvm_opt .= ", bed => q($OPTIONS{b})" if defined $OPTIONS{b};
				print PBS "JSNVM$idx=`qsub_command -m n -N ${sample}_${idx}_JSNVM -l walltime=$wt$depend -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::JointSNVMix; JointSNVMix(q($normal), q($tumor), [q(${vcf_out}_somatic_jsnvm.vcf)$germline], q($ref), { $jsnvm_opt });\\x27'`\n";
				$tumor{$sample}{vcf_comb}{jointsnvmix}{file} = "${vcf_out}_somatic_jsnvm.vcf";
				$tumor{$sample}{vcf_comb}{jointsnvmix}{depend} = ":".insert_version("JSNVM", "${vcf_out}_somatic_jsnvm.vcf", $sample, "\$JSNVM$idx", $idx);
			} elsif($caller eq "mutect") {
				# MuTect
				my $mut_opt = "Log => q($tumor{$sample}{Log}/${sample}_${idx}_Mut.txt)";
				$mut_opt .= ", germline_vcf => q(${vcf_out}_germline_mutect.vcf)" if $OPTIONS{g};
				$mut_opt .= ", cosmic => q($cosmic), dbsnp => q($dbsnp)" if $species eq "human";
				$mut_opt .= ", dbsnp => q($mouse_dbsnp)" if $species eq "mouse";
				$mut_opt .= ", dbsnp => q($fly_dbsnp)" if $species eq "fly";
				$mut_opt .= ", L => q($OPTIONS{b})" if defined $OPTIONS{b};
				print PBS "MUT$idx=`qsub_command -m n -N ${sample}_${idx}_Mut -l walltime=$wt$depend -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::GATK; muTect({normal => q($normal), tumor => q($tumor), somatic_vcf => q(${vcf_out}_somatic_mutect.vcf), R => q($ref), $mut_opt});\\x27'`\n";
				$tumor{$sample}{vcf_comb}{mutect}{file} = "${vcf_out}_somatic_mutect.vcf";
				$tumor{$sample}{vcf_comb}{mutect}{depend} = ":".insert_version("MUT", "${vcf_out}_somatic_mutect.vcf", $sample, "\$MUT$idx", $idx);
			} elsif($caller eq "id") {
				# IndelGenotyper
				my $id_opt = "Log => q($tumor{$sample}{Log}/${sample}_${idx}_ID.txt)";
				$id_opt .= ", germline_vcf => q(${vcf_out}_germline_id.vcf)" if $OPTIONS{g};
				$id_opt .= ", minFraction => 0.1, minCoverage => 6, minNormalCoverage => 4, minConsensusFraction => 0.7, minCnt => 1, window_size => 300";
				$id_opt .= ", L => q($OPTIONS{b})" if defined $OPTIONS{b};
				print PBS "ID$idx=`qsub_command -m n -N ${sample}_${idx}_ID -l walltime=$wt$depend -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::GATK; IndelGenotyper({normal => q($normal), tumor => q($tumor), somatic => q(${vcf_out}_indelgenotyper.vcf), R => q($ref), $id_opt});\\x27'`\n";
				$tumor{$sample}{vcf_comb}{indelgenotyper}{file} = "${vcf_out}_indelgenotyper.vcf";
				$tumor{$sample}{vcf_comb}{indelgenotyper}{depend} = ":".insert_version("ID", "${vcf_out}_indelgenotyper.vcf", $sample, "\$ID$idx", $idx);
			} elsif($caller eq "vs") {
				# Varscan
				$germline = ", q(${vcf_out}_germline_vs.vcf)" if $OPTIONS{g};
				my $vs_opt = "Log => q($tumor{$sample}{Log}/${sample}_${idx}_VS.txt)";
				$vs_opt .= ", bed => q($OPTIONS{b})" if defined $OPTIONS{b};
				print PBS "VS$idx=`qsub_command -m n -N ${sample}_${idx}_VS -l walltime=$wt$depend -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::Varscan; bam2vs_somatic2variant(q($normal), q($tumor), [q(${vcf_out}_somatic_loh_vs.vcf)$germline], q($ref), { $vs_opt });\\x27'`\n";

				if($OPTIONS{h}) {
					# remove loh calls into seperate file
					print PBS "VS_loh$idx=`qsub_command -m n -N ${sample}_${idx}_VS_LOH -W depend=afterok:\$VS$idx -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::Utils; filter_tsv_file(q(${vcf_out}_somatic_loh_vs.vcf), q(${vcf_out}_somatic_vs.vcf), { filter => \"[8] !~ \\\\\$;SS=3;\\\\\$\" });\\x27 > $tumor{$sample}{Log}/${sample}_${idx}_VS_LOH.txt 2>&1'`\n";
					print PBS "VS_Filt$idx=`qsub_command -m n -N ${sample}_${idx}_VS_Filt -W depend=afterok:\$VS$idx -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::Utils; filter_tsv_file(q(${vcf_out}_somatic_loh_vs.vcf), q(${vcf_out}_loh_vs.vcf), { filter => \"[8] =~ \\\\\$;SS=3;\\\\\$\" });\\x27 > $tumor{$sample}{Log}/${sample}_${idx}_VS_Filt.txt 2>&1'`\n";
					print PBS "RMVS$idx=`qsub_command -m n -N ${sample}_${idx}_VS_RM -W depend=afterok:\$VS_Filt$idx -j oe -d $OPTIONS{p}$queue 'rm -f ${vcf_out}_somatic_loh_vs.vcf'`\n";
					$tumor{$sample}{vcf_comb}{varscan}{file} = "${vcf_out}_somatic_vs.vcf";
					$tumor{$sample}{vcf_comb}{varscan}{depend} = ":".insert_version("VS_Filt", "${vcf_out}_somatic_vs.vcf", $sample, "\$VS_Filt$idx", $idx);
				} else {
					$tumor{$sample}{vcf_comb}{varscan}{file} = "${vcf_out}_somatic_loh_vs.vcf";
					$tumor{$sample}{vcf_comb}{varscan}{depend} = ":".insert_version("VS", "${vcf_out}_somatic_loh_vs.vcf", $sample, "\$VS$idx", $idx);
				}
			}
		}
		print PBS "\n";
		$idx++;
	}
	# combine vcfs
	combine_vcf(\%tumor);

	# add gatk variants if specified
	$idx = 1;
	$wt = int($walltime * 2 * 3600);
	foreach my $sample (sort keys %tumor) {
		if(defined $tumor{$sample}{vcf}{gatk} && @callers > 1) {
			my $vcf_out = $tumor{$sample}{vcf}{comb}{file};
			$vcf_out =~ s/\.vcf$/_ug.vcf/;
			print PBS "########### Checking to see if snps called in GATK for sample $sample\n";
			print PBS "SOM_GATK_Int$idx=`qsub_command -m n -N ${sample}_${idx}_SOM_GATK_Int -l walltime=$wt -W depend=afterok$tumor{$sample}{vcf}{comb}{depend}$tumor{$sample}{vcf}{gatk}{depend} -j oe -d $OPTIONS{p}$queue 'Somatic_GATK_intersect.pl -s $tumor{$sample}{vcf}{comb}{file} -o $vcf_out -g $tumor{$sample}{vcf}{gatk}{file} > $tumor{$sample}{Log}/${sample}_${idx}_SOM_GATK_Int.txt 2>&1'`\n";
			print PBS "RMSOM_GATK_Int$idx=`qsub_command -m n -N ${sample}_${idx}_RMSOM_GATK_Int -W depend=afterok:\$SOM_GATK_Int$idx -j oe -d $OPTIONS{p}$queue 'rm -f $tumor{$sample}{vcf}{comb}{file}*'`\n";
			$tumor{$sample}{vcf}{comb}{file} = $vcf_out;
			$tumor{$sample}{vcf}{comb}{depend} = ":\$RMSOM_GATK_Int$idx";
			$idx++;
		}
		print PBS "\n" if defined $tumor{$sample}{vcf}{gatk} && @callers > 1;
	}
}

# annotate
$idx = 1;
foreach my $sample (sort keys %tumor) {
	my $depend_all = "";
	foreach my $vcf (sort keys %{ $tumor{$sample}{vcf} }) {
		my $pbs_opt = "-j oe -d $OPTIONS{p}$queue";
		if($vcf eq 'hap' && $OPTIONS{c}) {
			next;
		} elsif($genome eq "human(hg19)" || $genome eq "mouse(mm9)" || $genome eq "mouse(mm10)" || $genome eq "fly(dm3)") {
			print PBS "########### Adding bam stats and Annotation for sample $sample and Vcf $vcf\n";
		} else {
			print PBS "########### Adding bam stats for sample $sample and Vcf $vcf\n";
			$pbs_opt .= " -M $email -m ae";
		}

		# add bam stats
		my $bastatsstring = join("), q(", @{ $tumor{$sample}{bamstats} });
		my $depend = $tumor{$sample}{vcf}{$vcf}{depend};
		my $wt = int($walltime * 2 * 3600);
		$depend .= $tumor{$sample}{depend} if defined $tumor{$sample}{depend};
		print PBS "BAM_${vcf}$idx=`qsub_command -N ${sample}_${idx}_${vcf}_BamStat -l walltime=$wt -W depend=afterok$depend $pbs_opt 'perl -e \\x27use Pipeline::VariantCalling::Annotation; add_bam_stats_to_vcf(q($ref), q($tumor{$sample}{vcf}{$vcf}{file}), [ q($bastatsstring) ], {Log => q($tumor{$sample}{Log}/${sample}_${idx}_${vcf}_BamStat.txt) });\\x27'`\n";
		$tumor{$sample}{vcf}{$vcf}{depend} = ":\$BAM_${vcf}$idx";

		if($genome eq "human(hg19)" || $genome eq "mouse(mm9)" || $genome eq "mouse(mm10)" || $genome eq "fly(dm3)") {
			my $vcf_base = fileparse($tumor{$sample}{vcf}{$vcf}{file}, ".vcf");
			my $info_rem = "$tumor{$sample}{TMP}/${vcf_base}_info_rem.vcf";
			my $vep = "$tumor{$sample}{VEP}/${vcf_base}_vep.tsv";
			my $tsv = "$tumor{$sample}{Var}/${vcf_base}.tsv";
	
			# remove the info column because sometimes this stuffs up vep (e.g. when using Pindel)
			$wt = int($walltime * 10 * 60);
			print PBS "Info_Rem_${vcf}$idx=`qsub_command -m n -N ${sample}_${idx}_${vcf}_Info_Rem -l walltime=$wt -W depend=afterok:\$BAM_${vcf}$idx $pbs_opt 'awk -v OFS=\"\\t\" \\x27{if(\$1 !~ /^#/) \$8 = \".\"; print }\\x27 $tumor{$sample}{vcf}{$vcf}{file} > $info_rem'`\n";

			# plugins
			my $plugins = "\"Additional_Annotation\", \"ExAC,/data/databases/ExAC/ExAC.r0.3_noTCGA.vcf.gz\", \"CADD,/data/databases/CADD/whole_genome_SNVs.tsv.gz,/data/databases/CADD/InDels.tsv.gz\"";
			$plugins .= ", \"Conservation\"" if $db_version == 73;
			$plugins .= ", \"EVS\"" if $species eq "human";
			if(defined $OPTIONS{e}) {
				$OPTIONS{e} =~ s/:/), q(/g;
				$plugins .= ", q($OPTIONS{e})";
			}
			print "Mouse(mm9) is only supported for genome build 67, overriding!!" if $genome eq "mouse(mm9)" && $db_version != 67;
			$db_version = 67 if $genome eq "mouse(mm9)";

			# run VEP
			$wt = int($walltime * 48 * 3600);
			print PBS "VEP_${vcf}$idx=`qsub_command -l nodes=1:ppn=4 -m n -N ${sample}_${idx}_${vcf}_VEP -l walltime=$wt -W depend=afterok:\$Info_Rem_${vcf}$idx $pbs_opt 'perl -e \\x27use Pipeline::VariantCalling::Annotation; variant_effect_predictor({everything => 1, input_file => q($info_rem), fork => 6, format => q(vcf), species => q($species), output_file => q($vep), check_existing => 1, check_alleles => 1, force_overwrite => 1, condel => \"b\", db_version => q($db_version), no_progress => 1, plugin => [ $plugins ]});\\x27 > $tumor{$sample}{Log}/${sample}_${idx}_${vcf}_Vep.txt 2>&1'`\n";
			#my $depend = insert_version("VEP", $vep, $sample, "\$VEP_${vcf}$idx", $idx);

			# fix html output
			$wt = int($walltime * 5 * 60);
			print PBS "HTML_${vcf}$idx=`qsub_command -m n -N ${sample}_${idx}_${vcf}_HTML -W depend=afterok:\$VEP_${vcf}$idx -l walltime=$wt $pbs_opt 'sed -i \"s/\\(chr_[^\\.]*\\)\\.\\([0-9a-zA-Z]*_area\\)/\\1_\\2/g\" ${vep}_summary.html'`\n";

			# convert to tsv
			my $env='';
			$env = "export PERL5LIB=\"/usr/local/cluster/all_arch/ensembl_api/api-ver-$db_version/ensembl-variation/modules\":\$PERL5LIB; " if defined $OPTIONS{o};
			my $TSV_options = "Log => q($tumor{$sample}{Log}/${sample}_${idx}_${vcf}_TSV.txt), db_version => $db_version";
			$TSV_options .= ", Vcf_quality => 1, Vcf_filter => 1" if $vcf eq 'comb';
			my $mail_opt = $OPTIONS{s}? '-m n' : "-M $email -m ae";
			$wt = int($walltime * 3600);
			print PBS "TSV_$vcf$idx=`qsub_command -N ${sample}_${idx}_${vcf}_TSV $mail_opt -l walltime=$wt -W depend=afterok:\$VEP_${vcf}$idx $pbs_opt '${env}perl -e \\x27use Pipeline::VariantCalling::Annotation; combine_vcf_vep(q($tumor{$sample}{vcf}{$vcf}{file}), q($vep), q($tsv), { $TSV_options });\\x27'`\n";
			$tumor{$sample}{vcf}{$vcf}{depend} = ":\$TSV_${vcf}$idx";

			if($OPTIONS{s}) {
				$wt = int($walltime * 3600);
				print PBS "Sort_$vcf$idx=`qsub_command -N ${sample}_${idx}_${vcf}_Sort -M $email -m ae -l walltime=$wt -W depend=afterok:\$TSV_${vcf}$idx $pbs_opt 'Sort_Variant_File.pl -s $sample \"$tsv\"'`\n";
				$tumor{$sample}{vcf}{$vcf}{depend} = ":\$Sort_${vcf}$idx";
			}
		} else {
			rmdir $tumor{$sample}{Var} if -d $tumor{$sample}{Var};
		}
		$depend_all .= $tumor{$sample}{vcf}{$vcf}{depend};
		print PBS "\n";
	}
	# remove tmp folder
	print PBS "RMVTMP$idx=`qsub_command -m n -N ${sample}_${idx}_RMVTMP -W depend=afterok$depend_all -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use File::Path qw(remove_tree); remove_tree q($tumor{$sample}{TMP});\\x27'`\n";
	print PBS "\n";
	$idx++;
}
print "Genome must be human(hg19), mouse(mm9), mouse(mm10) or fly(dm3) to work with the installed versions of the Ensembl Variant Effects Predictor, skipping.\n" unless $genome eq "human(hg19)" || $genome eq "mouse(mm9)" || $genome eq "mouse(mm10)" || $genome eq "fly(dm3)";
close(PBS);
system("sh $OPTIONS{p}/pbs_submit.sh") unless defined $OPTIONS{x};

sub ug {
	my $hash = shift;
	my $analysis_type = shift;
	my $idx = 1;

	print PBS "########### Running GATK Unified Genotyper\n";
	foreach my $sample (sort keys %{ $hash }) {
		# calculate vcf name
		my $vcf_out = fileparse($$hash{$sample}{current_bam_file});
		$vcf_out =~ s/\.bam$/_gatk.vcf/;
		$vcf_out =~ s/somatic_//;
		$vcf_out = "$$hash{$sample}{VCF}/$vcf_out";

		# pbs options
		my $wt = int($walltime * 24 * 3600);
		my $pbs_opt = "-l walltime=$wt -l nodes=1:ppn=8";
		$pbs_opt .= " -W depend=afterok$$hash{$sample}{depend}" if $$hash{$sample}{depend};

		# gatk options
		my $GATK_opt = "R => q($ref), nt => 4, nct => 2, minIndelFrac => 0";
		$GATK_opt .= ", dbsnp => q($dbsnp)" if $species eq "human";
		$GATK_opt .= ", dbsnp => q($mouse_dbsnp)" if $species eq "mouse";
		$GATK_opt .= ", dbsnp => q($fly_dbsnp)" if $species eq "fly";
		$GATK_opt .= ", L => q($OPTIONS{b})" if defined $OPTIONS{b};
		$GATK_opt .= ", Log => q($$hash{$sample}{Log}/${sample}_${idx}_UG.txt)";
	
		print PBS "UG$idx=`qsub_command -m n -N ${sample}_${idx}_UG $pbs_opt -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::GATK; UnifiedGenotyper({I => q($$hash{$sample}{current_bam_file}), o => q($vcf_out), $GATK_opt});\\x27'`\n";
		my $hashkey = "vcf_comb";
		$hashkey = 'vcf' if $analysis_type eq "tumor_norm";
		$$hash{$sample}{$hashkey}{gatk}{file} = $vcf_out;
		$$hash{$sample}{$hashkey}{gatk}{depend} = ":".insert_version("UG", $vcf_out, $sample, "\$UG$idx", $idx);
		$$hash{$sample}{depend_rm} = ":\$UG$idx";
		$idx++;
	}
	print PBS "\n";
}

sub pindel {
	my $hash = shift;
	my $idx = 1;

	print PBS "########### Running Pindel\n";
	foreach my $sample (sort keys %{ $hash }) {
		# calculate vcf name
		my $vcf_out = fileparse($$hash{$sample}{current_bam_file});
		$vcf_out =~ s/\.bam$/_pindel.vcf/;
		$vcf_out =~ s/somatic_//;
		$vcf_out = "$$hash{$sample}{VCF}/$vcf_out";

		# pbs options
		my $wt = int($walltime * 96 * 3600);
		my $pbs_opt = "-l walltime=$wt -l nodes=1:ppn=2";
		$pbs_opt .= " -W depend=afterok$$hash{$sample}{depend}" if $$hash{$sample}{depend};

		# pindel options
		my $pindel_opt = "fasta => q($ref), Log => q($$hash{$sample}{Log}/${sample}_${idx}_pindel.txt), number_of_threads => 2, TMP => q($$hash{$sample}{TMP})";
		$pindel_opt .= ", report_inversions => q(false), report_duplications => q(false), max_range_index => 1";
		$pindel_opt .= ", include => q($OPTIONS{b})" if defined $OPTIONS{b};

		print PBS "Pindel$idx=`qsub_command -m n -N ${sample}_${idx}_Pindel $pbs_opt -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::Pindel; Pindel([q(".join("), q(", @{ $tumor{$sample}{bamstats} }).")], q($vcf_out), {$pindel_opt});\\x27'`\n";

		# remove deletions with length > 100 base pairs as these won't get through VEP
		print PBS "Pindel_Rem_Lg_Del$idx=`qsub_command -m n -N ${sample}_${idx}_Pindel_Rem_Lg_Del -j oe -d $OPTIONS{p}$queue -W depend=afterok:\$Pindel$idx -l walltime=15:00 'sed -i \"/^[^#][^\\t]*\\t[0-9]*\\t\.\\t[ACTGN]\\{100\\}/d\" $vcf_out; igvtools index $vcf_out 2>&1'`\n";
		
		$$hash{$sample}{vcf}{pindel}{file} = $vcf_out;
		$$hash{$sample}{vcf}{pindel}{depend} = ":".insert_version("PINDEL", $vcf_out, $sample, "\$Pindel_Rem_Lg_Del$idx", $idx);
		$$hash{$sample}{depend_rm} = ":\$Pindel$idx";
		$idx++;
	}
	print PBS "\n";
}

sub hap {
	my $hash = shift;
	my $idx = 1;

	print PBS "########### Running GATK HaplotypeCaller\n";
	foreach my $sample (sort keys %{ $hash }) {
		# calculate vcf name
		my $vcf_out = fileparse($$hash{$sample}{current_bam_file});
		$vcf_out =~ s/\.bam$/_hap.vcf/;
		$vcf_out = "$$hash{$sample}{VCF}/$vcf_out";

		# pbs options
		my $wt = int($walltime * 72 * 3600);
		my $pbs_opt = "-l walltime=$wt -l nodes=1:ppn=8";
		$pbs_opt .= " -W depend=afterok$$hash{$sample}{depend}" if $$hash{$sample}{depend};

		# gatk options
		my $GATK_opt = "R => q($ref), nct => 8";
		$GATK_opt .= ", dbsnp => q($dbsnp)" if $species eq "human";
		$GATK_opt .= ", dbsnp => q($mouse_dbsnp)" if $species eq "mouse";
		$GATK_opt .= ", dbsnp => q($fly_dbsnp)" if $species eq "fly";
		$GATK_opt .= ", L => q($OPTIONS{b})" if defined $OPTIONS{b};
		$GATK_opt .= ", Log => q($$hash{$sample}{Log}/${sample}_${idx}_HAP.txt)";
		$GATK_opt .= ", emitRefConfidence => q(GVCF), variant_index_type => q(LINEAR), variant_index_parameter => 128000" if $OPTIONS{c};

		print PBS "HAP$idx=`qsub_command -m n -N ${sample}_${idx}_HAP $pbs_opt -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::GATK; HaplotypeCaller({I => q($$hash{$sample}{current_bam_file}), o => q($vcf_out), $GATK_opt});\\x27'`\n";
		my $hashkey = $OPTIONS{c}? "vcf":"vcf_comb";
		my $mail = $OPTIONS{c}? "-M $email -m ae ":"-m n ";
		$$hash{$sample}{$hashkey}{hap}{file} = $vcf_out;
		$$hash{$sample}{$hashkey}{hap}{depend} = ":".insert_version("HAP", $vcf_out, $sample, "\$HAP$idx", $idx, $mail);
		$idx++;
	}
	print PBS "\n";
}

sub vs_germ {
	my $hash = shift;
	my $idx = 1;

	print PBS "########### Running Varscan\n";
	foreach my $sample (sort keys %{ $hash }) {
		# calculate vcf name
		my $vcf_out = fileparse($$hash{$sample}{current_bam_file});
		$vcf_out =~ s/\.bam$/_vs.vcf/;
		$vcf_out = "$$hash{$sample}{VCF}/$vcf_out";

		# pbs options
		my $wt = int($walltime * 24 * 3600);
		my $pbs_opt = "-l walltime=$wt";
		$pbs_opt .= " -W depend=afterok$$hash{$sample}{depend}" if $$hash{$sample}{depend};

		# varscan options
		my $vs_opt = "Log => q($$hash{$sample}{Log}/${sample}_${idx}_VS.txt)";
		$vs_opt .= ", bed => q($OPTIONS{b})" if defined $OPTIONS{b};

		print PBS "VS$idx=`qsub_command -m n -N ${sample}_${idx}_VS $pbs_opt -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::Varscan; bam2mpileup2variant(q($$hash{$sample}{current_bam_file}), q($vcf_out), q($ref), { $vs_opt });\\x27'`\n";
		$$hash{$sample}{vcf_comb}{varscan}{file} = $vcf_out;
		$$hash{$sample}{vcf_comb}{varscan}{depend} = ":".insert_version("VS", $vcf_out, $sample, "\$VS$idx", $idx);
		$idx++;
	}
	print PBS "\n";
}

sub combine_vcf {
	my $hash = shift;
	my $idx = 1;

	foreach my $sample (sort keys %{ $hash }) {
		$$hash{$sample}{depend_rm} = "";

		if(keys %{ $$hash{$sample}{vcf_comb} } > 1) {
			print PBS "########### Combining vcfs $sample\n";
			my $depend;
			my @vcf_inputs;
			my @priority;

			my $vcf_out = fileparse($$hash{$sample}{current_bam_file});
			$vcf_out =~ s/\.bam$/_comb.vcf/;
			$vcf_out = "$$hash{$sample}{VCF}/$vcf_out";
			my $wt = int($walltime * 30 * 60);

			foreach my $vcf (sort keys %{ $$hash{$sample}{vcf_comb} }) {
				my $TMP_vcf = $$hash{$sample}{vcf_comb}{$vcf}{file};
				$TMP_vcf =~ s/\/VCF\//\/TMP\//;
				print PBS "UNIQ${idx}_$vcf=`qsub_command -m n -N ${sample}_${idx}_UNIQ_$vcf -l walltime=$wt -W depend=afterok$$hash{$sample}{vcf_comb}{$vcf}{depend} -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::Annotation; uniquify_vcf(q($$hash{$sample}{vcf_comb}{$vcf}{file}), q($vcf), {o => q($TMP_vcf)});\\x27'`\n";

				$depend .= ":\$UNIQ${idx}_$vcf";
				push @vcf_inputs, "q(:$vcf,vcf $TMP_vcf)";
				push @priority, $vcf;
			}

			$wt = int($walltime * 6 * 3600);
			print PBS "COMB$idx=`qsub_command -m n -N ${sample}_${idx}_COMB -l walltime=$wt -W depend=afterok$depend -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::GATK; CombineVariants({variant => [".join(", ", @vcf_inputs);
			print PBS "], minimalVCF => 1, R => q($ref), o => q($vcf_out), priority => q(".join(",", @priority )."), genotypemergeoption => \"PRIORITIZE\", setKey => \"Identified\", Log => q($$hash{$sample}{Log}/${sample}_${idx}_COMB.txt)});\\x27'`\n";
			$$hash{$sample}{vcf}{comb}{file} = $vcf_out;
			$$hash{$sample}{vcf}{comb}{depend} = ":\$COMB$idx";
			$$hash{$sample}{depend_rm} = ":\$COMB$idx";

			$idx++;
		} elsif(defined $$hash{$sample}{vcf}) {
			$$hash{$sample}{vcf_comb}{comb} = delete $$hash{$sample}{vcf_comb}{$_} foreach keys %{ $$hash{$sample}{vcf_comb} };
			my %new = ( %{ $$hash{$sample}{vcf} }, %{ $$hash{$sample}{vcf_comb} });
			%{ $$hash{$sample}{vcf} } = %new;
			$$hash{$sample}{depend_rm} .= $$hash{$sample}{vcf}{$_}{depend} foreach keys %{ $$hash{$sample}{vcf} };
		} else {
			%{ $$hash{$sample}{vcf} } = %{ $$hash{$sample}{vcf_comb} };
			$$hash{$sample}{depend_rm} .= $$hash{$sample}{vcf}{$_}{depend} foreach keys %{ $$hash{$sample}{vcf} };
		}
		print PBS "\n" if $idx > 1;
	}
}

sub recal {
	my $hash = shift;
	my $idx = 1;

	foreach my $sample (sort keys %{ $hash }) {
		print PBS "########### Performing GATK base quality recalibration on sample $sample\n";
		my $wt = int($walltime * 48 * 3600);
		my $pbs_opt = "-l walltime=$wt -l nodes=1:ppn=8";
		$pbs_opt .= " -W depend=afterok$$hash{$sample}{depend}" if defined $$hash{$sample}{depend};

		my $GATK_opt = "R => q($ref), nt => 8";
		$GATK_opt .= ", knownSites => q($dbsnp)" if $species eq "human";
		$GATK_opt .= ", knownSites => q($mouse_dbsnp)" if $species eq "mouse";
		$GATK_opt .= ", knownSites => q($fly_dbsnp)" if $species eq "fly";
		$GATK_opt .= ", L => q($OPTIONS{b})" if defined $OPTIONS{b};
		$GATK_opt .= ", Log => q($$hash{$sample}{Log}/${sample}_${idx}_Recal.txt)";
		my $newbam = $$hash{$sample}{current_bam_file};
		$newbam =~ s/\.bam$//;
		$newbam .= "_recal.bam";

		print PBS "RECAL$idx=`qsub_command -m n -N ${sample}_${idx}_Recal $pbs_opt -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::GATK; BaseRecalibrator({I => q($$hash{$sample}{current_bam_file}), o => q($newbam), $GATK_opt});\\x27'`\n";
		print PBS "RM$idx=`qsub_command -m n -N ${sample}_${idx}_Recal_RM -W depend=afterok:\$RECAL$idx -j oe -d $OPTIONS{p}$queue 'rm -f $$hash{$sample}{current_bam_file}'`\n" if $OPTIONS{i} || $OPTIONS{d} || $analysis_type ne "standard";

		$$hash{$sample}{current_bam_file} = $newbam;
		$$hash{$sample}{depend} = ":\$RECAL$idx";
		$idx++;
		print PBS "\n";
	}
}

sub indel {
	my $hash = shift;
	my $idx = 1;

	foreach my $sample (sort keys %{ $hash }) {
		print PBS "########### Running local indel aligner on sample $sample\n";
		my $wt = int($walltime * 48 * 3600);
		my $pbs_opt = "-l walltime=$wt -l nodes=1:ppn=8";
		$pbs_opt .= " -W depend=afterok$$hash{$sample}{depend}" if defined $$hash{$sample}{depend};

		my $GATK_opt = "R => q($ref), nt => 8";
		$GATK_opt .= ", known => [ q($mills), q($G1000_indels) ]" if $species eq "human";
		$GATK_opt .= ", L => q($OPTIONS{b})" if defined $OPTIONS{b};
		$GATK_opt .= ", Log => q($$hash{$sample}{Log}/${sample}_${idx}_Realign.txt)";
		my $newbam = $$hash{$sample}{current_bam_file};
		$newbam =~ s/\.bam$//;
		$newbam .= "_realign.bam";

		if($$BedFiles{$type}{aligner} eq "amplicon" || $$BedFiles{$type}{aligner} eq "amplicon_halo") {
			print PBS "INDELREALIGN$idx=`qsub_command -m n -N ${sample}_${idx}_Realign $pbs_opt -j oe -d $OPTIONS{p}$queue 'srma I=$$hash{$sample}{current_bam_file} O=$newbam R=$ref CREATE_INDEX=true MINIMUM_ALLELE_COVERAGE=1 MAXIMUM_TOTAL_COVERAGE=1000000 MAX_HEAP_SIZE=16384 > $$hash{$sample}{Log}/${sample}_${idx}_Realign.txt 2>&1'`\n";
		} else {
			print PBS "INDELREALIGN$idx=`qsub_command -m n -N ${sample}_${idx}_Realign $pbs_opt -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::VariantCalling::GATK; LocalIndelRealigner({I => q($$hash{$sample}{current_bam_file}), o => q($newbam), $GATK_opt});\\x27'`\n";
		}
		print PBS "RM$idx=`qsub_command -m n -N ${sample}_${idx}_Realign_RM -W depend=afterok:\$INDELREALIGN$idx -j oe -d $OPTIONS{p}$queue 'rm -f $$hash{$sample}{current_bam_file}'`\n" if $OPTIONS{d} || $analysis_type ne "standard";

		$$hash{$sample}{current_bam_file} = $newbam;
		$$hash{$sample}{depend} = ":\$INDELREALIGN$idx";
		$idx++;
		print PBS "\n";
	}
}

sub merge_bams {
	my $idx = shift;
	my $depends = shift;
	my $bams = shift;
	my $out = shift;
	my $log = shift;
	my $pbs_name = shift;
	my $input = "q(".join("), q(", @{ $bams }).")";
	my $depend = ${ $depends }[0];
	my $wt = int($walltime * 12 * 3600);
	my $pbs_opt = "-l walltime=$wt -l nodes=1:ppn=2 -l mem=12gb";
	$pbs_opt .= " -W depend=afterok".join("", @{ $depends }) if $depend;

	# merge the bam files
	print PBS "MERGE$idx=`qsub_command -m n -N $pbs_name $pbs_opt -j oe -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::Picard; MergeSamFiles([ $input ], q($out), { Log => q($log) });\\x27'`\n";
}

sub markdups {
	my $hash = shift;
	my $key = shift;
	my $name = shift;
	my $typean = $type;
	$typean =~ s/ /_/g if defined $typean;
	$typean = "_$typean" if defined $typean;
	my $idx = 1;
	my $wt = int($walltime * 12 * 3600);
	print PBS "########### Marking Duplicate Reads\n";
	foreach my $sample (sort keys %{ $hash }) {
		my $newbam = $$hash{$sample}{$key};
		$newbam =~ s/\.bam$//;
		$newbam .= "_dups.bam";
		print PBS "DUPS$name$idx=`qsub_command -m n -l walltime=$wt -j oe -N ${sample}_${idx}_Dups$name -d $OPTIONS{p}$queue 'perl -e \\x27use Pipeline::Picard; MarkDuplicates(q($$hash{$sample}{$key}), q($newbam), { Log => q($$hash{$sample}{Log}/${sample}_${idx}_Mark$typean.txt) });\\x27'`\n";
		$$hash{$sample}{depend} .= ":\$DUPS$name$idx";
		$$hash{$sample}{$key} = $newbam;
		$idx++;
	}
	print PBS "\n";
	return($idx);
}

sub symbolic_link {
	my $hash = shift;
	my $key = shift;
	print PBS "# Creating Symbolic Links\n";
	foreach my $sample(sort keys %{ $hash }) {
		# create a link to the bam file and add it to the bam stats command
		my ($bam, $path) = fileparse($$hash{$sample}{$key}, ".bam");
		$bam = "$$hash{$sample}{TMP}/$bam";
		my $bai = $$hash{$sample}{$key};
		$bai =~ s/m$/i/;
		print PBS "if [ ! -f \"$bam.bam\" ]; then\n\tln -s $$hash{$sample}{$key} $bam.bam\nfi\n";
		print PBS "if [ ! -f \"$bam.bai\" ]; then\n\tln -s $bai $bam.bai\nfi\n";
		$$hash{$sample}{"current_$key"} = "$bam.bam";
	}
	print PBS "\n";
}

sub insert_version {
	my $stage = shift;
	my $file = shift;
	my $sample = shift;
	my $depend = shift;
	my $idx = shift;
	my $mail = shift;
	$mail = "-m n " unless defined $mail;
	my $vers_info = "\#\#version_$Script=\\\"$version\\\"";
	$vers_info .= "\\n\#\#version_Picard_Mark_Duplicates=\\\"$PicardVerion\\\"" if $OPTIONS{d};
	$vers_info .= "\\n\#\#version_GATK_Indel_Realigner=\\\"$GATKversion\\\"" if $OPTIONS{i};
	$vers_info .= "\\n\#\#version_GATK_Base_Score_Recalibration=\\\"$GATKversion\\\"" if $OPTIONS{r};
	my $wt = int($walltime * 30 * 60);

	print PBS "VERS_$stage$idx=`qsub_command $mail -l walltime=$wt -j oe -N ${sample}_${idx}_Vers$stage -W depend=afterok:$depend -d $OPTIONS{p}$queue 'sed -i \"2i$vers_info\" $file'`\n";
	return("\$VERS_$stage$idx")
}

##code for manuscript "Evolution of late-stage metastatic melanoma is dominated by tetraploidization and aneuploidy"##

1) Variant Calling
2) Data organization
3) Post-processing

Code for Steps 1) and 2) provides execution steps for variant callers and folder structure generation for the output. These are included as a reference only for the variant calling done in Vergara et al, "Evolution of late-stage metastatic melanoma is dominated by tetraploidization and aneuploidy" as, formats, structures, etc will vary from person to person.   

Code for Step 3) includes post-processing steps, given the variant calling done in 1). Examples on how to execute each command is provided. 

1) Variant Calling
------------------

1.1) Run_Variant_Calling.pl

Given pre-processed bam files, runs Varscan, mutect and indelgenotyper.

1.2) create-run-mpileup.pl

Given pre-processed bam files, generates runs for mpileup.

1.3) create-multiSNV-chrom-level.pl
 
Given pre-processed bam files, generates per-chromosome runs for multisnv.

2) Data organization
--------------------

2.1) create_folder_structure.pl 

Generates the initial structure of folders that subsequent scripts will work on.

2.2) copy_variant_calling_output.pl 

Moves output from different variant callers into new folder structure.

2.3) build-multisnv-matrix-from-vcf.pl

Merges and reformats multisnv calls.

2.4) create-sample-vcf-for-multiSNV.pl  

Generates per-sample vcf for multisnv.

3) Post-processing
------------------ 

Steps provided with example output for patient MI-F.

PatientID: MI-F
Gender: Male (M)
Coverage Threshold: 10
Platform: WES
Germline ID: SRR2159460

3.3) generate-initial-merged-matrix.pl

Merges outputs from multisnv, mutect, indelgenotyper, varscan.

perl generate-initial-merged-matrix.pl MI-F SRR2159460 M 10 WES > MI-F_initial_merged_matrix.txt

3.4) filter-repeats-artifacts.R

Filters repetitive and artifactual regions in folder ./external_db

Rscript filter-repeats-artifacts.R  MI-F_initial_merged_matrix.txt  MI-F_initial_merged_matrix_filtReps.txt 

3.5) filter-non-overlapping-original-BED.R

Filters off-target reads

Rscript filter-non-overlapping-original-BED.R MI-F_initial_merged_matrix_filtReps.txt S03723424_Regions.bed MI-F_initial_merged_matrix_filtReps_withinBED.txt

3.6) generate-matrix-coverage.pl

Generates coverage matrix to be used to filter variants based on low coverage in any of the samples 

perl generate-matrix-coverage.pl MI-F_initial_merged_matrix_filtReps_withinBED.txt SRR2159460

3.7) filter-coverage.pl

Filteres based on coverage

perl filter-coverage.pl MI-F_initial_merged_matrix_filtReps_withinBED.txt MI-F_coverage_matrix.txt 10 > MI-F_initial_merged_matrix_filtReps_withinBED_10X.txt 

3.8) filter-GMAF-EVS.pl

Filters based on GMAF and EVS, resulting in the final matrix.

perl filter-GMAF-EVS.pl MI-F_initial_merged_matrix_filtReps_withinBED_10X.txt MI-F > MI-F.final.matrix.txt


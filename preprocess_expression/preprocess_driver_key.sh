################################################################################################################################################################
################################################################################################################################################################
#INPUT FILES
################################################################################################################################################################
#################################################################################################################################################################

#  Directory containing fastq files from 1st round of sequencing
fastq_round_1_input_dir="/project2/gilad/reem/heartproject/heart/fastq/"
#  Directory containing fastq files from 2nd round of sequencing
fastq_round_2_input_dir="/project2/gilad/reem/heartproject/heart2/fastq/"
#  Directory containing fastq files from sequencing of 4 additional samples (round 3)
fastq_round_3_input_dir="/project2/gilad/reem/heartproject/heart_4newlines/fastq/"
#  Directory containing fastq files from resequencing of 4 additional samples (round 4)
fastq_round_4_input_dir="/project2/gilad/reem/heartproject/heart_4newlines_reseq/fastq/"
# Directory containing fastq files from from 5 additional samples (round 5)
fastq_round_5_input_dir="/project2/gilad/reem/heartproject/heart_finalbatch/fastq/"
# Directory containing fastq files from those same 5 additional samples (round 6)
fastq_round_6_input_dir="/project2/gilad/reem/heartproject/heart_finalbatch_reseq/fastq/"


# The following three lane design files can be found here: https://github.com/relorbany/heart-pilot
# File used to map sequencing core names to (cell line, time step) names for the first round of sequencing. Produced by Katie Rhodes and Reem Elorbany
lane_design_round_1_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/LaneDesign2forR.csv"  
# File used to map sequencing core names to (cell line, time step) names for the second round of sequencing. Produced by Katie Rhodes and Reem Elorbany
lane_design_round_2_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/LaneDesignResequencing2forR.csv"
# File used to map sequencing core names to (cell line, time step) names for the third round of sequencing. Produced by Katie Rhodes and Reem Elorbany
# This is also used as reference for fastq_round_4
lane_design_round_3_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/LaneDesign_NewSamples.csv"

#  Genotype file created by Bryce Van de Geijn. It's format is vcf-ish.
#  Downloaded from "http://eqtl.uchicago.edu/jointLCL/" under the link "genotypes of 120 YRI individuals" on September 13, 2017
genotype_input="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/genotypesYRI.gen.txt.gz"

#  Heterozygous probabilities (from impute 2) genotype information contained in the following directory
#  Directory contains:
#       1. A file containing an ordered list of samples in the genotype files (YRI_samples.txt) 
#       2. One file for each chromosome of the file format "chr1.hg19.impute2.gz" that contains heterozygous probabilities for all snps in that chromosome
#  Bryce Van de Geijn sent me these files
heterozygous_site_input_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/heterozygous_probability_genotypes/"

#  Metadata/covariates compiled by Katie Rhodes and Reem Elorbany
# Can be downloaded from https://github.com/BennyStrobes/ipsc_cardiomyocyte_differentiation/tree/master/data/metadataUPDATED_04_11_2018.csv
metadata_input_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/metadataUPDATED_04_11_2018.csv"

#  Gencode hg19 gene annotation file
#  Downloaded from "https://www.gencodegenes.org/releases/19.html" on September 13, 2017
gencode_gene_annotation_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/gencode.v19.annotation.gtf.gz"

#  Directory containing files created by Nick Banovich and Bryce Van de Geijn
#  Data initially found on PPS cluster at /data/internal/genotypes/hg19/YRI/
impute2_genotype_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/impute2_genotypes/"

# Downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz on 10/20/17
# Required by WASP
chrom_info_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/chromInfo.txt"

# Directory containing some pre-compiled WASP functions
# Specifically, there are two important functions in this directory:
### 1. 'fasta2h5'
### 2. 'snp2h5'
# See https://github.com/bmvdgeijn/WASP/tree/master/snp2h5 for compiling these two functions
snp2h5_dir="/home/bstrober/ipsc_differentiation/preprocess/WASP/snp2h5/"

# Directory where tabix is installed
# Downloaded from https://sourceforge.net/projects/samtools/files/tabix/
tabix_directory="/project2/gilad/bstrober/tools/tabix-0.2.6/"

# Directory containing two files:
### 1. 'mixsvgp_K2_L100_1_28542829_assignments'
### 2. 'mixsvgp_K2_L20_0_29526888_gene_assignments'
### 3. 'mixsvgp_K2_L20_0_29526888_fmean'
### 4. 'flow_results.txt'
mixture_hmm_cell_line_grouping_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/mixture_hmm_cell_line_groupings/"

# File containing iPSC read counts from Banovich et al
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107654
ipsc_banovich_read_counts_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/expression_data_sets/iPSC_raw_counts.txt"

# File containing iPSC-CM read counts from Banovich et al
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107654
ipsc_cm_banovich_read_counts_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/expression_data_sets/GSE107654_iPSC-CM_counts.txt"


################################################################################################################################################################
################################################################################################################################################################
#OUTPUT DIRECTORIES (The following scripts assume these directories all exist)
################################################################################################################################################################
################################################################################################################################################################

# Output root for `preprocess_driver_key.sh`. All following directories are subdirectories of this.
preprocess_data_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/"

#  Directory containing merged (across sequencing rounds) fastq files for each sample
fastq_input_dir=$preprocess_data_dir"fastq/"

#  Directory containing fastqc results. One output file for every fastq input file
fastqc_dir=$preprocess_data_dir"fastqc_data/"

#  Directory containing reference genome (GRCh37)
genome_dir=$preprocess_data_dir"genome/"

#  Directory containing bams. One bam file for each sample.
bam_dir=$preprocess_data_dir"bam/"

#  Directory containing processed counts, quantile normalized expression data
preprocess_total_expression_dir=$preprocess_data_dir"processed_total_expression/"

#  Directory containing covariate information (covariates, PCs)
covariate_dir=$preprocess_data_dir"covariates/"

# Directory containing total expression data comparison to banovich et al
banovich_ipsc_comparison_dir=$preprocess_data_dir"banovich_ipsc_comparison/"

#  Directory containing plots/figures related to exploratory analysis of the total expression data (preprocess_total_expression_dir)
visualize_total_expression_dir=$preprocess_data_dir"visualize_total_expression/"

#  Directory containing various changes to $genotype_input so that it is ammendable to the software pipelines used to process allelic counts (ie WASP and GATK ASEReader)
genotype_dir=$preprocess_data_dir"genotype/"

#  Directory to contain various intermediate files developed by the WASP pipeline (many derivatives of the initial bams..)
wasp_intermediate_dir=$preprocess_data_dir"wasp_intermediate_files/"

#  Directory containing raw allelic counts
raw_allelic_counts_dir=$preprocess_data_dir"raw_allelic_counts/"

#  Directory containing processed allelic counts
processed_allelic_counts_dir=$preprocess_data_dir"processed_allelic_counts/"

#  Directory containing plots/figures related to exploratory analysis of the total expression data (preprocess_total_expression_dir)
visualize_allelic_counts_dir=$preprocess_data_dir"visualize_allelic_counts/"










#############################################################################################################################
# Scripts to preprocess total expresssion data
#############################################################################################################################

# Run the following 4 steps in series.
##########################################
##PART 1
# Concatenate all fastq files for each sample into one "merged" fastq file per sample. Samples have more than one fastq file initially because there are multiple sequencing rounds (to increase read depth)
# Also, create file called "$fastq_dir/fastq_mapping.txt" that contains mapping from merged fastq file to all of the corresponding replicate fastq file names. First column is merged fastq file name, second column is ','-seperated list of original files.
# Both merged fastq files and fastq_mapping.txt will be saved in output directory, $fastq_dir
# Note: This code really isn't the best b/c it very manually parses the fastq files based on the lane_design files. So caution should be taken in applying merge_fastq_replicates to new situations.
# Takes about 20 minutes to run
if false; then
sbatch merge_fastq_replicates.sh $fastq_round_1_input_dir $fastq_round_2_input_dir $fastq_round_3_input_dir $fastq_round_4_input_dir $fastq_round_5_input_dir $fastq_round_6_input_dir $lane_design_round_1_file $lane_design_round_2_file $lane_design_round_3_file $fastq_input_dir
fi


##Part 2
# Run "sbatch fastqc_and_download_reference_genome.sh $fastq_input_dir $fastqc_dir $genome_dir" to completion.
# This script runs fastqc on each of the fastq files, as well as downloads the reference genome
# Takes about 8 hours to run
if false; then
sbatch fastqc_and_download_reference_genome.sh $fastq_input_dir $fastqc_dir $genome_dir
fi


##PART 3
# Run "sh submit-subread.sh $fastq_input_dir $bam_dir $genome_dir". This shell submits (sbatch) a job for each sample (fastq file).
# Each job aligns fastq files using subread and creates and bam in $bam_dir
# Each job takes under 2 hours to run
if false; then
sh submit-subread.sh $fastq_input_dir $bam_dir $genome_dir
fi


# PART 4
# Run sbatch preprocess_total_expression.sh $preprocess_total_expression_dir $lane_design_file $exon_file $bam_dir $visualize_total_expression_dir $fastqc_dir $mixutre_hmm_cell_line_grouping_dir
# This script:
#    1. Processes the aligned read counts and creates quantile normalized expression data (preprocess_total_expression.R). 
#       # 1. also includes filters for genes. Genes need to have at least 10 samples such that RPKM >= .1 and raw counts >= 6
#    2. Re-order expression data into a matrix of dimension num_cell_lines X (genes*time_steps). This will be used for cell line PCA analysis (preprocess_total_expression_by_cell_lines.py)
#    3. Prepares covariate files (prepare_covariate_files.R)
#    4. Also does some exploratory visualization analysis of the expression data  (visualize_processed_total_expression.R)
#  Takes about 4 hours to run
exon_file=$genome_dir"exons.saf"
sh preprocess_total_expression.sh $preprocess_total_expression_dir $exon_file $bam_dir $visualize_total_expression_dir $metadata_input_file $covariate_dir $fastqc_dir $mixture_hmm_cell_line_grouping_dir $ipsc_banovich_read_counts_file $ipsc_cm_banovich_read_counts_file $banovich_ipsc_comparison_dir



#############################################################################################################################
# Preprocess allelic counts
#############################################################################################################################

# Run the following three parts in series:
#########################################



### Part 1: wasp_maping_pipeline_part1.sh. This includes:
######## A. Create text based SNP-site files (1 for each chromosome). These files only contain genotype sites (not information) and are used for WASP.
######## B. Make pseudo-VCF into VALID VCF file (does this for both vcf genotype file and heterozygous site file)
######## C. Get Reference genome in correct format for ASE Mapping
######## D. Convert impute2 genotype information to H5 format using WASP's snp2h5 script
######## E. Convert fasta information to h5 format using WASP's fasta2h5 script
sample_names=$fastq_input_dir"fastq_mapping.txt"
if false; then
sbatch wasp_mapping_pipeline_part1.sh $genotype_input $heterozygous_site_input_dir $genotype_dir $genome_dir $sample_names $impute2_genotype_dir $chrom_info_file $snp2h5_dir $tabix_directory
fi

# Genotype file that is in VCF format.
# Made by wasp_mapping_piepline_part1.sh
vcf_file=$genotype_dir"YRI_genotype.vcf.gz"



### Part 2: wasp_mapping_pipeline_part2.sh (run in parallel for each sample). This includes (per sample):
######## A. Map the fastq files using your favorite mapper/options and filter for quality using a cutoff of your choice (SUBREAD)
######## B. Use find_intersecting_snps.py (WASP script) to identify reads that may have mapping bias
######## C. Map the filtered fastq file a second time using the same arguments as the previous mapping
######## D. Use filter_remapped_reads.py (WASP script) to filter out reads where one or more of the allelic versions of the reads fail to map back to the same location as the original read
######## E. Merge the bams we plan to use. Then sort and index
######## F. Filter Duplicate Reads using rmdup.py (WASP script)
######## G. Run GATK ASE read counter (using output from wasp mapping pipeline)
######## H. Convert bam file from this individual into h5 format using bam2h5_tables_update.py (WASP script)
if false; then
while read standard_id_fastq sequencer_id; do
    standard_id=${standard_id_fastq::${#standard_id_fastq}-9}
    echo $standard_id
    sbatch wasp_mapping_pipeline_part2.sh $standard_id $genotype_dir $fastq_input_dir $wasp_intermediate_dir $genome_dir $vcf_file $raw_allelic_counts_dir $chrom_info_file
done<$sample_names
fi


































































#############################################
##OLD (retired scripts) NO LONGER USED!
#############################################


### PART 3
# Now that we have run the WASP mapping pipeline and GATK ASEReadCounter, we now have one allelic count file per sample
# process_and_organize_allelic_counts.sh will:
######### 1. Merge all of the sample specific count files into one table
######### 2. Map the sites to protein coding genes and remove sites that don't lie on a protein-coding gene
######### 3. For various heterozygous probability thresholds, place NA for (sample,site) pairs that have het. prob less than specified threshold
######### 4. Apply various filters for sites based on number of samples that we have mapped read counts to a het. site (etc)
######### 5. Visualize the number of counts we get at these various filters
if false; then
sbatch process_and_organize_allelic_counts.sh $raw_allelic_counts_dir $processed_allelic_counts_dir $genotype_dir $preprocess_total_expression_dir $gencode_gene_annotation_file $visualize_allelic_counts_dir
fi




# debug_sample_swap_driver.sh checks to make sure that every RNA-seq sample (has the correct label) and is paired correctly with its corresponding genotype
# We will test this by looping through each rna-seq sample, and for each sample:
##### Compare fraction of heterozygous sites that show bi-allelic expression for every possible genotype. 
##### There should only be one genotype that results in a high fraction. Make sure this genotype corresponds to the rna-seq sample
##### Takes about 8 hours to run.
if false; then
sbatch debug_sample_swap_driver.sh $processed_allelic_counts_dir $sample_swap_check_dir $genotype_dir $sample_names

# File containing expression files
expression_file_cross_data_sets="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/expression_data_sets/expression_samples.txt"

target_region_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/time_step_independent_qtl_pipelines/wasp/target_regions/target_regions_cis_distance_50000_maf_cutoff_0.1_min_reads_100_min_as_reads_25_min_het_counts_5_merged.txt"


##### Part 4
# Make heatmap showing correlation between ipsc temporal samples as well as those from other data sets
sh expression_correlation_heatmap.sh $expression_file_cross_data_sets $target_region_file $visualize_total_expression_dir 
fi


######################################
# TIME STEP INDEPENDENT WASP ANALYSIS
######################################

# This scripts assumes you have run the ipsc_preproccess_pipeline first (https://github.com/BennyStrobes/ipsc_preprocess_pipeline)
# And have save the results here:
preprocess_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/"



##########################################################################
# Input Files
##########################################################################
# Gencode hg19 gene annotation file
gencode_gene_annotation_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/gencode.v19.annotation.gtf.gz"

# Genotype directory created in preprocess pipeline (contains genotype data in h5 format)
genotype_dir=$preprocess_dir"genotype/"

# Directory created in preprocess pipeline that contains allelic count data in h5 format
raw_allelic_counts_dir=$preprocess_dir"raw_allelic_counts/"

# Downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz on 10/20/17
# Required by WASP
chrom_info_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/chromInfo.txt"


# Dosage genotype for all cell lines in our analysis
dosage_genotype_file=$preprocess_dir"genotype/YRI_genotype.vcf"


# Processed total expression data
# This is the quantile normalized expression data
quantile_normalized_expression=$preprocess_dir"processed_total_expression/quantile_normalized.txt"

# CM eqtl file
# Downloaded here http://eqtl.uchicago.edu/yri_ipsc/eQTL_WASP_CM.txt
cm_eqtl_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/eqtl_data_sets/eQTL_WASP_CM_thresh_1.0.txt"

# ipsc eqtl file
# Downloaded here http://eqtl.uchicago.edu/yri_ipsc/iPSC-eQTL-summary.txt
ipsc_eqtl_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/eqtl_data_sets/ipsc_eqtl_all/ipsc_eqtl_all_associations.txt"


##########################################################################
# Output Directories (Assumes these directories exist before script starts)
##########################################################################
# Root directory for ipsc qtl pipelines
wasp_qtl_root="/project2/gilad/bstrober/ipsc_differentiation_19_lines/time_step_independent_qtl_pipelines/wasp/"

# Directory containing target regions for WASP Test
target_regions_dir=$wasp_qtl_root"target_regions/"

# Directory containing CHT input files
cht_input_file_dir=$wasp_qtl_root"cht_input_files/"

# Directory containing CHT output files 
cht_output_dir=$wasp_qtl_root"cht_output/"

# Directory containing matrix factorization results
matrix_factorization_dir=$wasp_qtl_root"matrix_factorization/"

# Directory containing visualizations/plots
cht_visualization_dir=$wasp_qtl_root"cht_visualization/"



##########################################################################
# Model Parameters
##########################################################################
# Max distance between gene TSS and variant
cis_distance=50000

# Test variant must MAF >= maf_cutoff in each time step
maf_cutoff=.1

#  Must be >= $min_read_count (summed across all cell lines in 1 time step) reads at heterozygous variant!
min_read_count=100

# Must be >= $min_as_read_count (summed across all cell lines in 1 time step) reads on less popular (min transformation) allele
min_as_read_count=25

# Must be $min_het_count heterozygous variant in exon of gene of interest
min_het_count=5

# Keep track of model parameters
parameter_string="cis_distance_"$cis_distance"_maf_cutoff_0"$maf_cutoff"_min_reads_"$min_read_count"_min_as_reads_"$min_as_read_count"_min_het_counts_"$min_het_count







##################################################################
### Prepare input files to Combined Haplotype Test (CHT)
##################################################################
# Run the following 5 steps in series
##################################################################




###################################################################
### Step 1: get_wasp_target_regions.sh (run in parallel for each time step)
### This script makes a list of variant gene pairs (to be tested with CHT) for each time step that pass the following filters:
######## A. Distance between variant and gene TSS is <= $cis_distance
######## B. Test variant must have MAF >= $maf_cutoff
######## C. The gene must have >= $min_het_count heterozygous variants (summed across cell lines)
######## D. The number of reads overlapping heterozygous variants (summed across cell lines) must be >= $min_read_count
######## E. The number of reads overlapping heterozygous variants that map to less popular allele (min transformation) (summed across cell lines) must be >= $min_as_read_count
### Takes about 7 hours to run per time step
### This code was initially get_target_regions.py (in WASP repo), but I made some minor edits to it
if false; then
for time_step in $(seq 0 15); do
    sbatch get_wasp_target_regions.sh $time_step $target_regions_dir $quantile_normalized_expression $genotype_dir $raw_allelic_counts_dir $chrom_info_file $dosage_genotype_file $cis_distance $maf_cutoff $min_read_count $min_as_read_count $min_het_count $gencode_gene_annotation_file
done
fi






###################################################################
### Step 2: merge_wasp_target_regions.sh
### This script makes a list of variant gene pairs (to be tested with CHT) that pass step 1's filters in ALL time steps
### Basically just take the variant-gene pairs that are found in ALL time steps
### Takes about 15 minutes to run
if false; then
sbatch merge_target_regions.sh $target_regions_dir $cis_distance $maf_cutoff $min_read_count $min_as_read_count $min_het_count
fi

###################################################################
### Step 3: wasp_haplotype_shell.sh (run in parrallel for each sample)
### This script makes a CHT input file for each sample
### This script calls extract_haplotype_read_counts.py (WASP script)
### We run this script twice per sample:
###### 1. Once to create input file for WASP (homozygous test variants have no allelic read counts)
###### 2. Once to create input file for EAGLE2 (homozygous test variant do have allelic read counts)
### Takes about 8 hours per sample
if false; then
for time_step in $(seq 0 15); do
    sample_file=$target_regions_dir"rna_seq_samples_"$time_step".txt"
    while read cell_line; do
        sbatch wasp_haplotype_shell.sh $cell_line $time_step $genotype_dir $raw_allelic_counts_dir $chrom_info_file $parameter_string $target_regions_dir $cht_input_file_dir
    done<$sample_file
done
fi

##################################################################
# Step 4: Update total read depth (directly from WASP)
# This script updates the expected read depths of the input files for the Combined Haplotype Test. GC content biases and peakiness o
# of the data (potentially caused by differences in ChIP efficiency 
# or the expression of a few highly expressed genes) can vary greatly 
# between individuals, sequencing libraries, and lanes. To correct 
# for this, for each individual, we fit quartic splines across regions 
#for GC content and combined (across all individuals) read depth.
# The fitted splines provide an expected read depth for each
# region for each individual. This expected read depth is used 
# as an input to the CHT.
if false; then
sbatch update_total_depth.sh $genotype_dir $cht_input_file_dir
fi


###################################################################
### Step 5: fit_overdispersion_parameter.sh (run in parrallel for each time step)
### This script consists of two parts:
#### A. fit_as_coefficients.py: Estimate overdispersion parameters for allele-specific test (beta binomial)
#### B. fit_bnb_coefficients.py: Estimate overdispersion parameters for association test (beta-negative binomial)
#### C. get_PCs.R: Extract PCs based on haplotype read count data (make matrix num_samples by num_PCs where num_PCs == num_samples. Do it using quantile normalized expression)
# Takes about 3 hours per time step
if false; then
for time_step in $(seq 0 15); do
    sbatch fit_overdispersion_parameters.sh $time_step $parameter_string $cht_input_file_dir $target_regions_dir
done
fi









##################################################################
### Run Combined Haplotype Test (CHT) / WASP
##################################################################

# Run the following 4 steps in series
##################################################################

###################################################################
### Step 1: submit_chrom_parallel_cht_test.sh (run in parallel for each time step, chromosome, and number of pcs)
### This script runs cht for both real and permuted data
### This will produce MANY parallel jobs (6*16*22=2112)
### Each job takes about 5 hours to run (chromosome 1 takes a lot more!)
### WOWZERS!!!, thats a whole lot of compute hours.
### Sure is.
### DATA CRUNCHING!!!! NUM NUM NUM NUM.
### Because of this future, users probably shouldn't submit all at once.
if false; then
for pc_num in $(seq 0 5); do
    for time_step in $(seq 0 15); do 
        for chrom_num in $(seq 1 22); do 
            sbatch submit_chrom_parallel_cht_test.sh $time_step $chrom_num $pc_num $parameter_string $cht_input_file_dir $cht_output_dir
        done
    done
done
fi

###################################################################
### Step 2: organize_wasp_cht_test_results.sh (run in parallel for each time step)
### This script (does the following for each of the 6 possible numbers of PCs):
###  A. Concatenates the 22 chromosome WASP output files into 1 file (for both real and permuted data)
###  B. Run eFDR correction on data using real and permutation runs of wasp to assess genome wide significance
if false; then
for time_step in $(seq 0 15); do 
    sbatch organize_wasp_cht_test_results.sh $parameter_string $cht_input_file_dir $cht_output_dir $target_regions_dir $dosage_genotype_file $gencode_gene_annotation_file $quantile_normalized_expression $time_step $cis_distance
done
fi


##################################################################
### Step 3: Run down-stream analysis on WASP-CHT test results
### The following script runs three separate down-stream analyses for a range of PCs (0-5):
##### 1. `get_best_variant_per_gene.py`: Identify the most signficicant variant per egene (and gene with one significant per time step eQTL in any time step). The most significant variant is selected as the one with the smallest geometric mean pvalue across the 16 time steps.
##### 2.  `organize_tests_across_studies.py`: Using the egenes selected in part 1, find the pvalues of those variant gene pairs in Nick Banovich's iPSC and iPSC-CM eqtl data sets
##### 3. 'run_matrix_factorization.py': Run spare non-negative matrix factorization on the matrix of summary statisics (num_eGenesXnum_time_steps) for a range of number of latent factors and sparsity parameters
### It then runs `cht_visualization.R`: Make visualizations of WASP eqtl results
fdr=".05"
sh run_downstream_analysis_on_wasp_results.sh $parameter_string $cht_output_dir $fdr $target_regions_dir $matrix_factorization_dir $cm_eqtl_file $ipsc_eqtl_file $cht_visualization_dir




































































































######################################
# Old Analysis (not used currently)
#######################################


# CM eqtl egene file
cm_egene_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/eqtl_data_sets/eQTL_WASP_CM_egenes_thresh_1e-05.txt"

# ipsc eqtl egene file
ipsc_egene_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/eqtl_data_sets/permutations.all.RNAseq_run.fixed.txt.gz.bh.txt"



# Directory containing WASP/CHT summary statistics enriched in other data sets
cht_enrichment_dir=$wasp_qtl_root"enrichment/"

# File containing necessary information to convert from rsid to snpid and vice-versa
snp_id_to_rs_id_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/eqtl_data_sets/GTEx_Analysis_2016-01-15_v7_WGS_652ind_VarID_Lookup_Table.txt"


# File where each row contains information on an eqtl data set we want to compare with
# File was manually curated by me
eqtl_data_set_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/eqtl_data_sets/eqtl_data_sets_temp.txt"

# File containing egene and all-eqtl files for gtex tissues we are interested in comparing with
used_gtex_tissues_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/eqtl_data_sets/used_gtex_tissues.txt"

# Result of running metasoft on v7 eqtls
# Downloaded from: https://www.gtexportal.org/home/datasets on December 19, 2017
mvalue_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/eqtl_data_sets/GTEx_Analysis_v7.metasoft.txt"

# Directory containing chromHMM results
# Each cell line has its own file with suffix $cell_line_identifier'_15_coreMarks_mnemonics.bed.gz'
chrom_hmm_input_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/chrom_hmm/"



if false; then
python run_gene_set_enrichment_analysis.py $parameter_string $cht_output_dir $cht_visualization_dir $target_regions_dir $gencode_gene_annotation_file $ipsc_egene_file $ipsc_eqtl_file $cm_egene_file $cm_eqtl_file $eqtl_data_set_file
fi


num_permutations="1000"
efdr=".1"
if false; then
for pc_num in $(seq 0 5); do
    if false; then
    for time_step in $(seq 0 15); do 
        sbatch cre_enrichment_analysis.sh $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $efdr $time_step
    done
    fi
    Rscript visualize_cre_enrichment_analysis.R $parameter_string $cht_visualization_dir $pc_num $cht_enrichment_dir $num_permutations $efdr 
done
fi

if false; then
sh banovich_cre_enrichment_analysis.sh $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $cm_egene_file $num_permutations $ipsc_eqtl_file $cm_eqtl_file $parameter_string $cht_enrichment_dir $cht_visualization_dir
fi





if false; then
sbatch cluster_variant_gene_pairs.sh $parameter_string $cht_output_dir $cht_visualization_dir
fi


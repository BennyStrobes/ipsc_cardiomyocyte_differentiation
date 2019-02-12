#!/bin/bash
#SBATCH --time=10:00:00 --partition=broadwl --mem=5GB

parameter_string="$1"
qtl_results_dir="$2"
real_eqtl_results_file="$3"
perm_eqtl_results_file="$4"
num_jobs="$5"


# Merge seperate dynamic eQTL runs into one file for real data
python merge_seperate_dynamic_eqtl_runs.py $real_eqtl_results_file $qtl_results_dir $parameter_string"_permute_False_results" $num_jobs
# Merge seperate dynamic eQTL runs into one file for permuted data
python merge_seperate_dynamic_eqtl_runs.py $perm_eqtl_results_file $qtl_results_dir $parameter_string"_permute_True_results" $num_jobs



################
# eFDR analysis
# output file for eFDR analysis
efdr_file=$qtl_results_dir$parameter_string"_eFDR_results.txt"
# Run eFDR correction
Rscript eFDR_correction.R $real_eqtl_results_file $perm_eqtl_results_file $efdr_file

################
# Assess genome wide significance of actual data based on the eFDR approach with FDR <= $fdr_thresh
# Output file for all significant variant gene pairs
fdr_thresh=".05"
significant_efdr_results=$qtl_results_dir$parameter_string"_efdr_"$fdr_thresh"_significant.txt"
# Output file for significant egenes and their strongest associated variant
significant_efdr_gene_results=$qtl_results_dir$parameter_string"_efdr_"$fdr_thresh"_significant_egenes.txt"
python assess_significance_efdr_approach.py $efdr_file $real_eqtl_results_file $significant_efdr_results $significant_efdr_gene_results $fdr_thresh

################
# Assess genome wide significance of actual data based on the eFDR approach with FDR <= $fdr_thresh
# Output file for all significant variant gene pairs
fdr_thresh=".01"
# Assess genome wide significance of actual data based on the eFDR approach with FDR <= $fdr_thresh
# Output file for all significant variant gene pairs
significant_efdr_results=$qtl_results_dir$parameter_string"_efdr_"$fdr_thresh"_significant.txt"
# Output file for significant egenes and their strongest associated variant
significant_efdr_gene_results=$qtl_results_dir$parameter_string"_efdr_"$fdr_thresh"_significant_egenes.txt"
python assess_significance_efdr_approach.py $efdr_file $real_eqtl_results_file $significant_efdr_results $significant_efdr_gene_results $fdr_thresh


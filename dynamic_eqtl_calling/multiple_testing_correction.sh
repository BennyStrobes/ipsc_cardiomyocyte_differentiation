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
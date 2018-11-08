#!/bin/bash
#SBATCH --time=10:00:00 --partition=broadwl --mem=5GB

model_version="$1"
covariate_method="$2"
num_jobs="$3"
parameter_string="$4"
dynamic_eqtl_input_file="$5"
qtl_results_dir="$6"


real_eqtl_results_file=$qtl_results_dir$parameter_string"_permute_False_results.txt"
perm_eqtl_results_file=$qtl_results_dir$parameter_string"_permute_True_results.txt"
sh multiple_testing_correction.sh $parameter_string $qtl_results_dir $real_eqtl_results_file $perm_eqtl_results_file $num_jobs
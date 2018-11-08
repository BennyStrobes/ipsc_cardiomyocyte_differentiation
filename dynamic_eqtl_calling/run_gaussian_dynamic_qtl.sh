#!/bin/bash
#SBATCH --time=10:00:00 --partition=broadwl --mem=5GB

input_data_file="$1"
output_file="$2"
model_version="$3"
permute="$4"
covariate_method="$5"
job_number="$6"
num_jobs="$7"


num_lines=`wc -l $input_data_file`

# Run dynamic qtl 
Rscript run_gaussian_dynamic_qtl.R $input_data_file $output_file $model_version $permute $covariate_method $job_number $num_jobs $num_lines



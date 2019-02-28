#!/bin/bash
#SBATCH --time=4:00:00 --mem=10GB --partition=broadwl
time_step="$1"
parameter_string="$2"
cht_input_file_dir="$3"
target_regions_dir="$4"

CHT_IN_FILE=$cht_input_file_dir"cht_input_file_"$parameter_string"_time_"$time_step".txt"
ls $cht_input_file_dir"haplotype_read_counts_"$parameter_string*"_"$time_step".txt.gz" | grep -v adjusted > $CHT_IN_FILE
date

#
# Estimate overdispersion parameters for allele-specific test (beta binomial)
#
AS_COEF_OUT_FILE=$cht_input_file_dir"cht_as_coef_"$parameter_string"_time_"$time_step".txt"
python fit_as_coefficients.py $CHT_IN_FILE $AS_COEF_OUT_FILE


date
#
# Estimate overdispersion parameters for association test (beta-negative binomial)
#
BNB_COEF_OUT_FILE=$cht_input_file_dir"cht_bnb_coef_"$parameter_string"_time_"$time_step".txt"
python fit_bnb_coefficients.py --min_counts 50 --min_as_counts 10 $CHT_IN_FILE $BNB_COEF_OUT_FILE

date




samples_file=$target_regions_dir"rna_seq_samples_"$time_step".txt"
PC_output_file=$cht_input_file_dir"pcs_"$parameter_string"_time_"$time_step".txt"
Rscript get_PCs.R $samples_file $cht_input_file_dir $parameter_string $time_step $PC_output_file

samples_file=$target_regions_dir"rna_seq_samples_"$time_step".txt"
PC_pve_output_file=$cht_input_file_dir"pcs_"$parameter_string"_time_"$time_step"_variance_explained.txt"
Rscript get_PCs_variance_explained.R $samples_file $cht_input_file_dir $parameter_string $time_step $PC_pve_output_file
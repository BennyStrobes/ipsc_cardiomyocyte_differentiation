#!/bin/bash
#SBATCH --time=1:00:00 --partition=broadwl --mem=5GB

joint_test_input_file="$1"
dynamic_eqtl_input_file="$2"
total_expression_file="$3"
genotype_file="$4"
environmental_variable_form="$5"
genotype_version="$6"
cell_line_specific_pc_file="$7"
target_region_input_file="$8"

python preprocess_data_for_dynamic_eQTL_calling.py $joint_test_input_file $dynamic_eqtl_input_file $total_expression_file $genotype_file $environmental_variable_form $genotype_version $cell_line_specific_pc_file $target_region_input_file

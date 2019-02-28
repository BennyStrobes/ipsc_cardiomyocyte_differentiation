#!/bin/bash
#SBATCH --time=1:00:00 --mem=10GB --partition=broadwl

parameter_string="$1"
cht_output_dir="$2"
fdr="$3"
target_regions_dir="$4"
matrix_factorization_dir="$5"
cm_eqtl_file="$6"
ipsc_eqtl_file="$7"
cht_visualization_dir="$8"
cht_input_file_dir="$9"

if false; then
for pc_num in $(seq 0 5); do

	echo $pc_num

	# Part 1: Identify the most signficicant variant per egene (and gene with one significant per time step eQTL in any time step). The most significant variant is selected as the one with the smallest geometric mean pvalue across the 16 time steps.
	python get_best_variant_per_gene.py $parameter_string $cht_output_dir $pc_num $fdr

	# Part 2: Using the egenes selected in part 1, find the pvalues of those variant gene pairs in Nick Banovich's iPSC and iPSC-CM eqtl data sets
	python organize_tests_across_studies.py $parameter_string $cht_output_dir $pc_num $cm_eqtl_file $ipsc_eqtl_file $fdr

	# Part 3: Run spare non-negative matrix factorization on the matrix of summary statisics (num_eGenesXnum_time_steps) for a range of number of latent factors and sparsity parameters
	python run_matrix_factorization.py $parameter_string $cht_output_dir $matrix_factorization_dir $target_regions_dir $fdr $pc_num

done
fi

# Make visualizations of WASP eqtl results
Rscript cht_visualization.R $parameter_string $cht_output_dir $cht_visualization_dir $matrix_factorization_dir $cht_input_file_dir
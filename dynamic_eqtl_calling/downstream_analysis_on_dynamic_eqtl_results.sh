#!/bin/bash
#SBATCH --time=10:00:00 --partition=broadwl --mem=5GB

model_version="$1"
covariate_method="$2"
num_jobs="$3"
parameter_string="$4"
dynamic_eqtl_input_file="$5"
qtl_results_dir="$6"
qtl_pvalue_distribution_visualization_dir="$7"
cell_line_overlap_analysis_dir="$8"
genotype_file="$9"
t_0_eqtl_results_file="${10}"


########################################
### Part A: Multiple testing correction
### Merges results from parallelization runs and computes significance after multiple testing correction
real_eqtl_results_file=$qtl_results_dir$parameter_string"_permute_False_results.txt"
perm_eqtl_results_file=$qtl_results_dir$parameter_string"_permute_True_results.txt"
if false; then
sh multiple_testing_correction.sh $parameter_string $qtl_results_dir $real_eqtl_results_file $perm_eqtl_results_file $num_jobs
fi

########################################
### Part B: Visualize dynamic eQTL pvalue distributions
### Plot pvalue distributions (qq-plots) for dynamic eqtl run
if false; then
Rscript visualize_dynamic_eqtl_pvalue_distribution.R $real_eqtl_results_file $perm_eqtl_results_file $parameter_string $qtl_pvalue_distribution_visualization_dir
fi



########################################
### Part C: Cell Line overlap analysis
# For each cell line pair, compute fraction of time (across dynamic eQTLs and background variants) that those two cell lines were in the same genotype bin ({0,1,2})
# Do so for both dynamic eQTLs and per-time step eqtls at time step 0
if false; then
num_genes="200"
real_overlap_matrix=$cell_line_overlap_analysis_dir$parameter_string"_"$num_genes"_genes_real_overlap_matrix.txt"
perm_overlap_matrix=$cell_line_overlap_analysis_dir$parameter_string"_"$num_genes"_genes_perm_overlap_matrix.txt"
python perform_cell_line_overlap_analysis.py $real_eqtl_results_file $genotype_file $real_overlap_matrix $perm_overlap_matrix $num_genes "dynamic_eqtl"
real_overlap_matrix=$cell_line_overlap_analysis_dir$parameter_string"_"$num_genes"_genes_t0_real_overlap_matrix.txt"
perm_overlap_matrix=$cell_line_overlap_analysis_dir$parameter_string"_"$num_genes"_genes_t0_perm_overlap_matrix.txt"
python perform_cell_line_overlap_analysis.py $t_0_eqtl_results_file $genotype_file $real_overlap_matrix $perm_overlap_matrix $num_genes "standard_eqtl"
fi
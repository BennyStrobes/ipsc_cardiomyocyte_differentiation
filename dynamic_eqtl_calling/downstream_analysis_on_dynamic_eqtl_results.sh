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
time_step_independent_stem="${10}"
chrom_hmm_input_dir="${11}"
tissue_specific_chrom_hmm_enrichment_dir="${12}"
time_step_independent_comparison_dir="${13}"
gsea_data_dir="${14}"
gencode_file="${15}"
gene_set_enrichment_dir="${16}"
cardiomyopathy_gene_list="${17}"
gtex_gwas_hits_dir="${18}"
gwas_overlap_dir="${19}"
liftover_directory="${20}"
visualization_input_dir="${21}"

date

########################################
### Part A: Multiple testing correction
### Merges results from parallelization runs and computes significance after multiple testing correction
real_eqtl_results_file=$qtl_results_dir$parameter_string"_permute_False_results.txt"
perm_eqtl_results_file=$qtl_results_dir$parameter_string"_permute_True_results.txt"
sh multiple_testing_correction.sh $parameter_string $qtl_results_dir $real_eqtl_results_file $perm_eqtl_results_file $num_jobs
significant_egene_file=$qtl_results_dir$parameter_string"_efdr_.05_significant_egenes.txt"
significant_qtl_file=$qtl_results_dir$parameter_string"_efdr_.05_significant.txt"


########################################
### Part B: Visualize dynamic eQTL pvalue distributions
### Plot pvalue distributions (qq-plots) for dynamic eqtl run
Rscript visualize_dynamic_eqtl_pvalue_distribution.R $real_eqtl_results_file $perm_eqtl_results_file $parameter_string $qtl_pvalue_distribution_visualization_dir




########################################
### Part C: Cell Line overlap analysis
# For each cell line pair, compute fraction of time (across dynamic eQTLs and background variants) that those two cell lines were in the same genotype bin ({0,1,2})
# Do so for both dynamic eQTLs and per-time step eqtls at time step 0
num_genes="200"
real_overlap_matrix=$cell_line_overlap_analysis_dir$parameter_string"_"$num_genes"_genes_real_overlap_matrix.txt"
perm_overlap_matrix=$cell_line_overlap_analysis_dir$parameter_string"_"$num_genes"_genes_perm_overlap_matrix.txt"
python perform_cell_line_overlap_analysis.py $real_eqtl_results_file $genotype_file $real_overlap_matrix $perm_overlap_matrix $num_genes "dynamic_eqtl"
real_overlap_matrix=$cell_line_overlap_analysis_dir$parameter_string"_"$num_genes"_genes_t0_real_overlap_matrix.txt"
perm_overlap_matrix=$cell_line_overlap_analysis_dir$parameter_string"_"$num_genes"_genes_t0_perm_overlap_matrix.txt"
python perform_cell_line_overlap_analysis.py $time_step_independent_stem"0_eqtl_results.txt" $genotype_file $real_overlap_matrix $perm_overlap_matrix $num_genes "standard_eqtl"


########################################
### Part D: Tissue specific chromHMM enrichment analysis
# Compute enrichment of dynamic eQTLs within cell type matched chromHMM enhancer elements
threshold="1.0"
num_permutations="1000"
sh tissue_specific_chrom_hmm_enrichment_analysis.sh $parameter_string $real_eqtl_results_file $significant_egene_file $num_permutations $threshold $chrom_hmm_input_dir $time_step_independent_stem $model_version $tissue_specific_chrom_hmm_enrichment_dir




########################################
### Part E: Time Step Independent Comparison
# Compare Dynamic eQTLs to per time step eQTLs
threshold="1.0"
dynamic_standard_egenes_comparison_file=$time_step_independent_comparison_dir$parameter_string"_"$threshold"_dynamic_standard_egenes_comparison.txt"
dynamic_standard_egenes_background_comparison_file=$time_step_independent_comparison_dir$parameter_string"_"$threshold"_background_dynamic_standard_egenes_comparison.txt"
python time_step_independent_comparison.py $dynamic_standard_egenes_comparison_file $dynamic_standard_egenes_background_comparison_file $time_step_independent_stem $threshold $model_version $real_eqtl_results_file $significant_egene_file




########################################
### Part F: Gene Set enrichment within GSEA
python gsea_gene_set_enrichment_analysis.py $parameter_string $significant_egene_file $gencode_file $gene_set_enrichment_dir $gsea_data_dir $time_step_independent_stem 


########################################
### Part G: Gene Set enrichment within dilated cardiomyopathy gene sets
python cardiomyopathy_gene_set_enrichment_analysis.py $parameter_string $significant_egene_file $gencode_file $gene_set_enrichment_dir $cardiomyopathy_gene_list $parameter_string $real_eqtl_results_file


########################################
### Part H: Enrichment within GTEx GWAS variants
threshold="5e-8"
python gtex_gwas_dynamic_qtl_overlap.py $gtex_gwas_hits_dir $gwas_overlap_dir $significant_qtl_file $parameter_string $threshold
threshold="5e-6"
python gtex_gwas_dynamic_qtl_overlap.py $gtex_gwas_hits_dir $gwas_overlap_dir $significant_qtl_file $parameter_string $threshold


########################################
### Part I: Extract GWAS data for Miami plots at a few specific, exemplary positions
python extract_specific_gwas_examples_for_miami_plots.py $gtex_gwas_hits_dir $gwas_overlap_dir$parameter_string"_" $real_eqtl_results_file $genotype_file $liftover_directory

########################################
### Part J: Organize significant eqtl results for dynamic eqtl visualization

python organize_dynamic_qtl_egenes_for_visualization.py $dynamic_eqtl_input_file $significant_qtl_file $visualization_input_dir$parameter_string

date


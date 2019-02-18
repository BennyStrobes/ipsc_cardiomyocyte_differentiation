#!/bin/bash
#SBATCH --time=10:00:00 --partition=broadwl --mem=5GB



qtl_results_dir="$1"
cell_line_overlap_analysis_dir="$2"
tissue_specific_chrom_hmm_enrichment_dir="$3"
time_step_independent_comparison_dir="$4"
gwas_overlap_dir="$5"
visualization_dir="$6"
dynamic_eqtl_input_file="$7"


Rscript visualize_dynamic_qtl_results.R $qtl_results_dir $cell_line_overlap_analysis_dir $tissue_specific_chrom_hmm_enrichment_dir $time_step_independent_comparison_dir $gwas_overlap_dir $visualization_dir $dynamic_eqtl_input_file
#!/bin/bash
#SBATCH --time=10:00:00 --partition=broadwl --mem=5GB


parameter_string="$1"
real_eqtl_results_file="$2"
significant_egene_file="$3"
num_permutations="$4"
threshold="$5"
chrom_hmm_input_dir="$6"
time_step_independent_stem="$7"
model_version="$8"
tissue_specific_chrom_hmm_enrichment_dir="$9"


hits_versions=( "early_time_step_hits" "late_time_step_hits" "change_in_sign_hits")

for hits_version in "${hits_versions[@]}"; do
	#######################
	marker_type="enhancer"
	#######################
	cell_line_version="heart_only_cell_lines"
	mod_parameter_string=$parameter_string"_"$marker_type"_"$cell_line_version"_"$hits_version"_"$num_permutations"_"$threshold
	python perform_tissue_specific_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $real_eqtl_results_file $significant_egene_file $time_step_independent_stem $tissue_specific_chrom_hmm_enrichment_dir$mod_parameter_string $hits_version $model_version $threshold

	cell_line_version="ipsc_only_cell_lines"
	mod_parameter_string=$parameter_string"_"$marker_type"_"$cell_line_version"_"$hits_version"_"$num_permutations"_"$threshold
	python perform_tissue_specific_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $real_eqtl_results_file $significant_egene_file $time_step_independent_stem $tissue_specific_chrom_hmm_enrichment_dir$mod_parameter_string $hits_version $model_version $threshold
done

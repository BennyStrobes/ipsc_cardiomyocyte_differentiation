#!/bin/bash
#SBATCH --time=18:00:00 --mem=10GB

time_step="$1"
target_regions_dir="$2"
corrected_quantile_normalized_expression="$3"
genotype_dir="$4"
raw_allelic_counts_dir="$5"
chrom_info_file="$6"
dosage_genotype_file="$7"
cis_distance="$8"
maf_cutoff="$9"
min_read_count="${10}"
min_as_read_count="${11}"
min_het_count="${12}"
gencode_gene_annotation_file="${13}"


date

# First need to make list of all the cell lines in this time step
cell_lines_in_current_time_step=$target_regions_dir"rna_seq_samples_"$time_step".txt"
python get_cell_line_names_in_each_time_step.py $time_step $corrected_quantile_normalized_expression $cell_lines_in_current_time_step


python get_wasp_target_regions.py $corrected_quantile_normalized_expression $gencode_gene_annotation_file $dosage_genotype_file $genotype_dir $cis_distance $maf_cutoff $raw_allelic_counts_dir $min_read_count $min_as_read_count $chrom_info_file $time_step $min_het_count $target_regions_dir

date



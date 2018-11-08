#!/bin/bash
#SBATCH --time=02:30:00 --mem=20GB --partition=broadwl

raw_allelic_counts_dir="$1"
processed_allelic_counts_dir="$2"
genotype_dir="$3"
preprocess_total_expression_dir="$4"
gencode_gene_annotation_file="$5"
visualize_allelic_counts_dir="$6"

# Now that we have run the WASP mapping pipeline and GATK ASEReadCounter, we now have one allelic count file per sample
# process_and_organize_allelic_counts.sh will:
######### 1. Merge all of the sample specific count files into one table
######### 2. Map the sites to protein coding genes and remove sites that don't lie on a protein-coding gene
######### 3. For various heterozygous probability thresholds, place NA for (sample,site) pairs that have het. prob less than specified threshold
######### 4. Apply various filters for sites based on number of samples that we have mapped read counts to a het. site (etc)
python process_and_organize_allelic_counts.py $raw_allelic_counts_dir $processed_allelic_counts_dir $genotype_dir $preprocess_total_expression_dir $gencode_gene_annotation_file



Rscript visualize_processed_allelic_counts.R $processed_allelic_counts_dir $visualize_allelic_counts_dir $genotype_dir $preprocess_total_expression_dir


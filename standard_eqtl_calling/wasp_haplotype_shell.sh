#!/bin/bash
#SBATCH --time=22:00:00 --mem=15GB --partition=broadwl

INDIVIDUAL="$1"
time_step="$2"
genotype_dir="$3"
raw_allelic_counts_dir="$4"
chrom_info_file="$5"
parameter_string="$6"
target_regions_dir="$7"
cht_input_file_dir="$8"



echo $INDIVIDUAL
echo $time_step

ALL_SAMPLES_FILE=$genotype_dir"all_genotyped_samples.txt"
target_region_file=$target_regions_dir"target_regions_"$parameter_string"_merged.txt"

date
#
# create CHT input file for this individual (For WASP)
# Create input file for WASP (homozygous test variants have no allelic read counts)
#
python extract_haplotype_read_counts.py \
    --chrom $chrom_info_file \
    --snp_index $genotype_dir"snp_index.h5" \
    --snp_tab $genotype_dir"snp_tab.h5" \
    --geno_prob $genotype_dir"geno_probs.h5" \
    --haplotype $genotype_dir"haps.h5" \
    --samples $ALL_SAMPLES_FILE \
    --individual $INDIVIDUAL \
    --ref_as_counts $raw_allelic_counts_dir"ref_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --alt_as_counts $raw_allelic_counts_dir"alt_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --other_as_counts $raw_allelic_counts_dir"other_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --read_counts $raw_allelic_counts_dir"read_counts."$INDIVIDUAL"_"$time_step".h5" \
    $target_region_file \
    | gzip > $cht_input_file_dir"haplotype_read_counts_"$parameter_string"."$INDIVIDUAL"_"$time_step".txt.gz"


#
# create CHT input file for this individual (For EAGLE 2)
# Create input file for EAGLE2 (homozygous test variants do have allelic read counts)
#
python extract_haplotype_read_counts.py \
    --chrom $chrom_info_file \
    --snp_index $genotype_dir"snp_index.h5" \
    --snp_tab $genotype_dir"snp_tab.h5" \
    --geno_prob $genotype_dir"geno_probs.h5" \
    --haplotype $genotype_dir"haps.h5" \
    --samples $ALL_SAMPLES_FILE \
    --individual $INDIVIDUAL \
    --ref_as_counts $raw_allelic_counts_dir"ref_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --alt_as_counts $raw_allelic_counts_dir"alt_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --other_as_counts $raw_allelic_counts_dir"other_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --read_counts $raw_allelic_counts_dir"read_counts."$INDIVIDUAL"_"$time_step".h5" \
    --homozygous_as_counts "rand_hap" \
    $target_region_file \
    | gzip > $cht_input_file_dir"haplotype_read_counts_rand_hap_"$parameter_string"."$INDIVIDUAL"_"$time_step".txt.gz"

date
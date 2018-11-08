#!/bin/bash
#SBATCH --time=1:00:00

target_regions_dir="$1"
cis_distance="$2"
maf_cutoff="$3"
min_read_count="$4"
min_as_read_count="$5"
min_het_count="$6"

date
python merge_target_regions.py $target_regions_dir $cis_distance $maf_cutoff $min_read_count $min_as_read_count $min_het_count
date
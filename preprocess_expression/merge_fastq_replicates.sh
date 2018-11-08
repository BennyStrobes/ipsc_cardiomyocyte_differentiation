#!/bin/bash
#SBATCH --time=01:00:00


fastq_round_1_input_dir="$1"
fastq_round_2_input_dir="$2"
fastq_round_3_input_dir="$3"
fastq_round_4_input_dir="$4"
fastq_round_5_input_dir="$5"
fastq_round_6_input_dir="$6"
lane_design_round_1_file="$7"
lane_design_round_2_file="$8"
lane_design_round_3_file="${9}"
fastq_dir="${10}"

python merge_fastq_replicates.py $fastq_round_1_input_dir $fastq_round_2_input_dir $fastq_round_3_input_dir $fastq_round_4_input_dir $fastq_round_5_input_dir $fastq_round_6_input_dir $lane_design_round_1_file $lane_design_round_2_file $lane_design_round_3_file $fastq_dir
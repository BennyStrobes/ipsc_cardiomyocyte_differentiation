#!/bin/bash
#SBATCH --time=26:00:00 --mem=15GB --partition=broadwl


genotype_dir="$1"
cht_input_file_dir="$2"


seq_file=$genotype_dir"seq.h5"


module unload Anaconda2

IN_FILE=$cht_input_file_dir"input_files.txt"
OUT_FILE=$cht_input_file_dir"output_files.txt"
ls $cht_input_file_dir"haplotype_read_counts_cis_"* | grep -v adjusted > $IN_FILE
cat $IN_FILE | sed 's/.txt/.adjusted.txt/' >  $OUT_FILE



read_depth_coefficient_file=$cht_input_file_dir"read_depth_coefficients.txt"
python update_total_depth.py --skip 40 --seq $seq_file --fit_out_file $read_depth_coefficient_file $IN_FILE $OUT_FILE



python update_total_depth.py --seq $seq_file --fit_in_file $read_depth_coefficient_file $IN_FILE $OUT_FILE


source ~/.bash_profile



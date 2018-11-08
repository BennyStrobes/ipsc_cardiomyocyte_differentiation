#!/bin/bash
#SBATCH --time=26:00:00 --partition=broadwl

time_step="$1"
chrom_num="$2"
pc_num="$3"
parameter_string="$4"
cht_input_file_dir="$5"
cht_output_pc_opti_dir="$6"

module unload Anaconda3
module load python


echo $time_step
echo $chrom_num
echo $pc_num

CHT_IN_FILE=$cht_input_file_dir"cht_input_file_"$parameter_string"_time_"$time_step".txt"
AS_COEF_OUT_FILE=$cht_input_file_dir"cht_as_coef_"$parameter_string"_time_"$time_step".txt"
BNB_COEF_OUT_FILE=$cht_input_file_dir"cht_bnb_coef_"$parameter_string"_time_"$time_step".txt"
PC_output_file=$cht_input_file_dir"pcs_"$parameter_string"_time_"$time_step".txt"


date

#
# run combined haplotype test on real data
#
CHT_OUT_FILE=$cht_output_pc_opti_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_"$chrom_num".txt"
python combined_test_parrallel_chromosome.py --min_as_counts 0 \
       --bnb_disp $BNB_COEF_OUT_FILE \
       --as_disp $AS_COEF_OUT_FILE \
       --num_pcs $pc_num --pc_file $PC_output_file \
       --chrom_num $chrom_num \
       $CHT_IN_FILE $CHT_OUT_FILE

date

#
# run combined haplotype test on permuted data
#
CHT_OUT_FILE=$cht_output_pc_opti_dir"cht_perm1_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_"$chrom_num".txt"
python combined_test_parrallel_chromosome.py --min_as_counts 0 \
       --bnb_disp $BNB_COEF_OUT_FILE \
       --as_disp $AS_COEF_OUT_FILE \
       --num_pcs $pc_num --pc_file $PC_output_file \
       --chrom_num $chrom_num \
       --shuffle \
       $CHT_IN_FILE $CHT_OUT_FILE

date

if false; then
CHT_OUT_FILE=$cht_output_pc_opti_dir"cht_perm2_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_"$chrom_num".txt"
python combined_test_parrallel_chromosome.py --min_as_counts 0 \
       --bnb_disp $BNB_COEF_OUT_FILE \
       --as_disp $AS_COEF_OUT_FILE \
       --num_pcs $pc_num --pc_file $PC_output_file \
       --chrom_num $chrom_num \
       --shuffle \
       $CHT_IN_FILE $CHT_OUT_FILE
fi
echo "DONE"

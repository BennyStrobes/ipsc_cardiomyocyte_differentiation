#!/bin/bash
#SBATCH --time=15:00:00 --mem=10GB --partition=broadwl


parameter_string="$1"
cht_input_file_dir="$2"
cht_output_dir="$3"
target_regions_dir="$4"
dosage_genotype_file="$5"
gencode_gene_annotation_file="$6"
corrected_quantile_normalized_expression="$7"
time_step="$8"
cis_distance="$9"





# Some parameters to this function
target_region_file=$target_regions_dir"target_regions_"$parameter_string"_merged.txt"


# Do these for varying number of PCs
for pc_num in $(seq 0 5); do

    ###################################################
    # Concatenate all 22 chromosome's results file into one output
    ###################################################
    echo $pc_num"_"$time_step
    # Concatenate all 22 chromosomes of real data
    # Also, modify/subset some of the columns from WASP output to make the output easier to handle
    wasp_results_stem=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_"
    python organize_wasp_output.py $time_step $cht_output_dir $dosage_genotype_file $corrected_quantile_normalized_expression $gencode_gene_annotation_file $target_region_file $wasp_results_stem
    # Concatenate all 22 chromosomes of permuted data
    # Also, modify/subset some of the columns from WASP output to make the output easier to handle    
    wasp_results_stem=$cht_output_dir"cht_perm1_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_"
    python organize_wasp_output.py $time_step $cht_output_dir $dosage_genotype_file $corrected_quantile_normalized_expression $gencode_gene_annotation_file $target_region_file $wasp_results_stem


    # Result files from concatenation
    null_file1=$cht_output_dir"cht_perm1_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eqtl_results.txt"
    # null_file2=$cht_output_dir"cht_perm2_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eqtl_results.txt"
    real_file=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eqtl_results.txt"
    
    ###################################################
    # Run eFDR correction
    ###################################################
    efdr_file=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eFDR_results.txt"
    Rscript eFDR_correction.R $real_file $null_file1 $efdr_file

    fdr_thresh=".05"
    # Assess genome wide significance of actual data based on the eFDR approach with FDR <= $fdr_thresh
    significant_efdr_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_efdr_thresh_"$fdr_thresh"_significant.txt"
    significant_efdr_gene_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_efdr_thresh_"$fdr_thresh"_significant_egenes.txt"
    python assess_wasp_significance_efdr_approach.py $efdr_file $real_file $significant_efdr_results $significant_efdr_gene_results $fdr_thresh


    fdr_thresh=".1"
    # Assess genome wide significance of actual data based on the eFDR approach with FDR <= $fdr_thresh
    significant_efdr_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_efdr_thresh_"$fdr_thresh"_significant.txt"
    significant_efdr_gene_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_efdr_thresh_"$fdr_thresh"_significant_egenes.txt"
    python assess_wasp_significance_efdr_approach.py $efdr_file $real_file $significant_efdr_results $significant_efdr_gene_results $fdr_thresh


    fdr_thresh=".2"
    # Assess genome wide significance of actual data based on the eFDR approach with FDR <= $fdr_thresh
    significant_efdr_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_efdr_thresh_"$fdr_thresh"_significant.txt"
    significant_efdr_gene_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_efdr_thresh_"$fdr_thresh"_significant_egenes.txt"
    python assess_wasp_significance_efdr_approach.py $efdr_file $real_file $significant_efdr_results $significant_efdr_gene_results $fdr_thresh




    ###############
    # OLD!
    ###############


    if false; then 



    # Using reference eqtl data sets, extract pvalues belonging to reference variant-gene pairs from OUR DATA
    eqtl_file=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eqtl_results.txt"
    enrichment_output_stem=$cht_enrichment_dir"enrichment_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step
    # python get_snp_gene_pairs_from_other_data_and_match_for_background.py $eqtl_file $enrichment_output_stem $eqtl_data_set_file $cis_distance $mvalue_file



    # Assess genome wide significance of actual data based on the emperical FDR threshold $fdr_thresh
    significant_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_fdr_thresh_"$fdr_thresh"_significant.txt"
    significant_gene_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_fdr_thresh_"$fdr_thresh"_significant_egenes.txt"
    #python assess_wasp_significance.py $null_file $real_file $significant_results $significant_gene_results $fdr_thresh


    # Compute qvalues on null data
    null_file=$cht_output_dir"cht_perm_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eqtl_results.txt"
    real_file=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eqtl_results.txt"
    qvalue_file=$cht_output_dir"null_qvalue_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step".txt"
    #Rscript compute_qvalues.R $null_file $qvalue_file

    # Assess genome wide significance of actual data based on the qvalue threshold (FDR <= .1) in null data
    significant_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_qval_"$qval_thresh"_significant.txt"
    significant_gene_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_qval_"$qval_thresh"_significant_egenes.txt"
    #python assess_wasp_significance.py $null_file $real_file $qvalue_file $significant_results $significant_gene_results $qval_thresh

    fi
done
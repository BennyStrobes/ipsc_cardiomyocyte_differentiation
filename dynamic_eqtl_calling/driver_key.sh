
######################################
# Dynamic eQTL analysis
######################################

# This scripts assumes you have run the ipsc_preproccess_pipeline first (https://github.com/BennyStrobes/ipsc_cardiomyocyte_differentiation/tree/master/preprocess_expression)
# And have save the results at the root directory here:
preprocess_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/"

# This scripts assumes you have run the ipsc_standard_eqtl_pipeline first (https://github.com/BennyStrobes/ipsc_cardiomyocyte_differentiation/tree/master/standard_eqtl_calling)
# And have save the results at the root directory here:
standard_eqtl_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/time_step_independent_qtl_pipelines/wasp/"





###############################################################################
# Input Data
###############################################################################

# File containing all of the target regions we are using. We are using this file to convert from gene positions to ensable id
target_region_input_file=$standard_eqtl_dir"target_regions/target_regions_cis_distance_50000_maf_cutoff_0.1_min_reads_100_min_as_reads_25_min_het_counts_5_merged.txt"

# cell line specific pcs
cell_line_specific_pc_file=$preprocess_dir"covariates/cell_line_ignore_missing_principal_components_10.txt"

# Gene expression data for all samples

total_expression_file=$preprocess_dir"processed_total_expression/quantile_normalized_no_projection.txt"

# Dosage-based genotypes for all samples
genotype_file=$preprocess_dir"genotype/YRI_genotype.vcf"

# Directory containing chromHMM results
# Each cell line has its own file with suffix $cell_line_identifier'_15_coreMarks_mnemonics.bed.gz'
chrom_hmm_input_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/chrom_hmm/"

# Time step independent data
# Suffixes are:
#### 1. "$time_step"_eqtl_results.txt, containing pvalues for each time step, along with MAF, and dist-toTSS for our variant-gene pairs
#### 2. "$time_step"_efdr_thresh_.1_significant_egenes.txt giving list of significant variant gene pairs at this time step
time_step_independent_stem="/project2/gilad/bstrober/ipsc_differentiation_19_lines/time_step_independent_qtl_pipelines/wasp/cht_output/cht_results_cis_distance_50000_maf_cutoff_0.1_min_reads_100_min_as_reads_25_min_het_counts_5_num_pc_3_time_"

# Directory containing gsea data (for gene set enrichment)
gsea_data_dir="/project2/gilad/bstrober/tools/tools/gsea/data/"

# File containing conversions from ensamble ids to gene symbols
gencode_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/gencode.v19.annotation.gtf.gz"

# File containing gwas hits
gwas_catalog="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/gwas_catalog_v1.0.2-associations_e92_r2018-07-17.tsv"

# Directory containing gtex gwas hits
gtex_gwas_hits_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/gtex_gwas_data/"

# File containing list of cardiomyopathy genes
cardiomyopathy_gene_list="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/cardiomyopathy_disease_genes/cardiomyopathy_genes.txt"

# Directory contaiing executables required to run liftover
liftover_directory="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/liftOver_x86/"

# iPSC-CM eqtl file created by Nick Banovich et al (https://genome.cshlp.org/content/early/2017/12/05/gr.224436.117)
# Downloaded here http://eqtl.uchicago.edu/yri_ipsc/eQTL_WASP_CM.txt
cm_eqtl_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/eqtl_data_sets/eQTL_WASP_CM_thresh_1.0.txt"

# iPSC eqtl file created by Nick Banovich et al (https://genome.cshlp.org/content/early/2017/12/05/gr.224436.117)
# Downloaded here http://eqtl.uchicago.edu/yri_ipsc/iPSC-eQTL-summary.txt
ipsc_eqtl_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/eqtl_data_sets/ipsc_eqtl_all/ipsc_eqtl_all_associations.txt"



###############################################################################
# Output directories (aasume all of these exist prior to starting analysis)
###############################################################################

# Root directory for this of all ipsc data based results
output_root="/project2/gilad/bstrober/ipsc_differentiation_19_lines/gaussian_dynamic_qtl_pipelines_v2/"

# Directory containing necessary input files to qtl tests
input_data_dir=$output_root"input_data/"

# Directory containing text files with results from dynamic qtl analysis
temp_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/gaussian_dynamic_qtl_pipelines/"
qtl_results_dir=$temp_dir"qtl_results/"

# Directory containing visualization of results found in qtl_results_dir
qtl_pvalue_distribution_visualization_dir=$output_root"qtl_pvalue_distribution_visualization/"

# Directory containing results from cell_line_overlap_analysis
cell_line_overlap_analysis_dir=$output_root"cell_line_overlap/"

# Output directory containing tissue specific chromHMM enrichment results
tissue_specific_chrom_hmm_enrichment_dir=$output_root"tissue_specific_chromHMM_enrichment/"

# Output directory containing comparisons to time-step independent analysis
time_step_independent_comparison_dir=$output_root"time_step_independent_comparison/"

# Output directory containing gene set enrichments
gene_set_enrichment_dir=$output_root"gene_set_enrichment/"

# Output directory containing gwas overlaps
gwas_overlap_dir=$output_root"gwas_overlap/"

# Output directory containing comparison other eqtl data sets
eqtl_data_set_comparison_dir=$output_root"eqtl_data_set_comparison/"

# Output directory containing input data for dynamic eqtl visualiztion
visualization_input_dir=$output_root"visualization_input/"

# Output directory containing manuscript figures
visualization_dir=$output_root"visualization/"

# Output directory containing results from power analysis
power_analysis_dir=$output_root"power_analysis/"






##########################################
# Step 1: Preprocess data for dynamic eQTL calling
##########################################
# prepare files necessary to run dynamic eQTL calling
# Two files necessary for dynamic eQTL calling (and are created in this `preprocess_data_for_dynamic_eQTL_calling.sh`:
###### 1. $joint_test_input_file: File with one line for each tested sample. Each line contains sample name, environmental variable, and first 8 cell line PCs for that sample
###### 2. $dynamic_eqtl_input_file: File with one line for each variant-gene pair tested in dynamic eQTL analysis. Each line contains lists of expression, genotype, cell line PCs, and environmental variables across all tested samples
# Note: Takes about 30 minutes to run

# How to encode environmental variable:
### 1. "time_steps": Encode environmental variable with samples day of differentiation
environmental_variable_form="time_steps"
# How to encode genotype data
### 1. "dosage": Use dosage based genotype values
### 2. "round": Round dosage based genotype values to {0,1,2]}
genotype_version="dosage"

# First output file fron this analysis
joint_test_input_file=$input_data_dir"joint_test_input_file_"$environmental_variable_form".txt"
# Second output file from this analysis
dynamic_eqtl_input_file=$input_data_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_input.txt"
if false; then
sbatch preprocess_data_for_dynamic_eQTL_calling.sh $joint_test_input_file $dynamic_eqtl_input_file $total_expression_file $genotype_file $environmental_variable_form $genotype_version $cell_line_specific_pc_file $target_region_input_file
fi









##########################################
# Step 2: Run GLM/GLMM dynamic qtl modeling
##########################################

##########################
# Parameters
##########################
## 1. $model_version: Takes on "glm", "glm_quadratic", or "glmm"
## 2. covariate_method: Takes on "none", "pc1", "pc1_2", "pc1_3", "pc1_4", "pc1_5"
## 3. num_jobs: How many nodes to parallelize on  to parrallelize

num_jobs="10"

################################
# Run Dynamic eQTLs using GLM with sweep over covariate methods and model versions
covariate_methods=( "none" "pc1" "pc1_2" "pc1_3" "pc1_4" "pc1_5")
model_versions=( "glm" "glmm" "glm_quadratic")
################################

# Loop through covariate methods
for covariate_method in "${covariate_methods[@]}"; do
    # Loop through model versions
    for model_version in "${model_versions[@]}"; do
        # Loop through parallelization jobs
        for job_number in $(seq 0 $(($num_jobs-1))); do 
            
            # Run for real data
            permute="False"
            output_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_"$permute"_results_"$job_number".txt"
            if false; then
            sbatch run_gaussian_dynamic_qtl.sh $dynamic_eqtl_input_file $output_file $model_version $permute $covariate_method $job_number $num_jobs
            fi
            # Run for permuted data
            permute="True"
            output_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_"$permute"_results_"$job_number".txt"
            if false; then
            sbatch run_gaussian_dynamic_qtl.sh $dynamic_eqtl_input_file $output_file $model_version $permute $covariate_method $job_number $num_jobs
            fi
        done
    done
done

# Run for a few more models
covariate_methods=( "pc1_6" "pc1_7" "pc1_8" "pc1_9" "pc1_10")
model_versions=( "glm" )
# Loop through covariate methods
for covariate_method in "${covariate_methods[@]}"; do
    # Loop through model versions
    for model_version in "${model_versions[@]}"; do
        # Loop through parallelization jobs
        for job_number in $(seq 0 $(($num_jobs-1))); do 
            
            # Run for real data
            permute="False"
            output_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_"$permute"_results_"$job_number".txt"
            if false; then
            sbatch run_gaussian_dynamic_qtl.sh $dynamic_eqtl_input_file $output_file $model_version $permute $covariate_method $job_number $num_jobs
            fi
            # Run for permuted data
            permute="True"
            output_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_"$permute"_results_"$job_number".txt"
            if false; then
            sbatch run_gaussian_dynamic_qtl.sh $dynamic_eqtl_input_file $output_file $model_version $permute $covariate_method $job_number $num_jobs
            fi
        done
    done
done




##########################################
# Step 3: Run Downstream analysis on eQTL results
##########################################

# The following script runs many types of downstream analysis on the dynamic eqtl results:
########################################
### Part A: Multiple testing correction
### Merges results from parallelization runs and computes significance after multiple testing correction
########################################
### Part B: Visualize dynamic eQTL pvalue distributions
### Plot pvalue distributions (qq-plots) for dynamic eqtl run
########################################
### Part C: Cell Line overlap analysis
### For each cell line pair, compute fraction of time (across dynamic eQTLs and background variants) that those two cell lines were in the same genotype bin ({0,1,2})
### Do so for both dynamic eQTLs and per-time step eqtls at time step 0
########################################
### Part D: Tissue specific chromHMM enrichment analysis
# Compute enrichment of dynamic eQTLs within cell type matched chromHMM enhancer elements
########################################
### Part E: Time Step Independent Comparison
# Compare Dynamic eQTLs to per time step eQTLs
########################################
### Part F: Gene Set enrichment within GSEA
########################################
### Part G: Gene Set enrichment within dilated cardiomyopathy gene sets
########################################
### Part H: Enrichment within GTEx GWAS variants
########################################
### Part I: Extract GWAS data for Miami plots at a few specific, exemplary positions
########################################
### Part J: Organize significant eqtl results for dynamic eqtl visualization
########################################
### Part K: Compare dynamic eqtls to existing data sets


covariate_methods=( "none" "pc1" "pc1_2" "pc1_3" "pc1_4" "pc1_5")
model_versions=( "glm" "glmm" "glm_quadratic")


# Loop through covariate methods
if false; then
for covariate_method in "${covariate_methods[@]}"; do
    # Loop through model versions
    for model_version in "${model_versions[@]}"; do
        parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
        sbatch downstream_analysis_on_dynamic_eqtl_results.sh $model_version $covariate_method $num_jobs $parameter_string $dynamic_eqtl_input_file $qtl_results_dir $qtl_pvalue_distribution_visualization_dir $cell_line_overlap_analysis_dir $genotype_file $time_step_independent_stem $chrom_hmm_input_dir $tissue_specific_chrom_hmm_enrichment_dir $time_step_independent_comparison_dir $gsea_data_dir $gencode_file $gene_set_enrichment_dir $cardiomyopathy_gene_list $gtex_gwas_hits_dir $gwas_overlap_dir $liftover_directory $visualization_input_dir $cm_eqtl_file $ipsc_eqtl_file $eqtl_data_set_comparison_dir
    done
done
fi

covariate_methods=( "pc1_6" "pc1_7" "pc1_8" "pc1_9" "pc1_10")
model_versions=( "glm" )
if false; then
for covariate_method in "${covariate_methods[@]}"; do
    # Loop through model versions
    for model_version in "${model_versions[@]}"; do
        parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
        sbatch downstream_analysis_on_dynamic_eqtl_results.sh $model_version $covariate_method $num_jobs $parameter_string $dynamic_eqtl_input_file $qtl_results_dir $qtl_pvalue_distribution_visualization_dir $cell_line_overlap_analysis_dir $genotype_file $time_step_independent_stem $chrom_hmm_input_dir $tissue_specific_chrom_hmm_enrichment_dir $time_step_independent_comparison_dir $gsea_data_dir $gencode_file $gene_set_enrichment_dir $cardiomyopathy_gene_list $gtex_gwas_hits_dir $gwas_overlap_dir $liftover_directory $visualization_input_dir $cm_eqtl_file $ipsc_eqtl_file $eqtl_data_set_comparison_dir
    done
done
fi




##########################################
# Step 4: Perform dynamic qtl power analysis
##########################################
if false; then
sh run_power_analysis.sh $power_analysis_dir
fi

##########################################
# Step 5: Visualize results from dynamic qtl analysis
##########################################
sh visualize_dynamic_qtl_results.sh $qtl_results_dir $cell_line_overlap_analysis_dir $tissue_specific_chrom_hmm_enrichment_dir $time_step_independent_comparison_dir $gwas_overlap_dir $eqtl_data_set_comparison_dir $visualization_input_dir $visualization_dir $power_analysis_dir













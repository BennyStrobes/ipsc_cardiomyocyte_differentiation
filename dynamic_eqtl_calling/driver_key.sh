
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
cell_line_specific_pc_file=$preprocess_dir"covariates/cell_line_ignore_missing_principal_components_9.txt"

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




# Output directory containing visualizations of top dynanmic qtl visualizations
dynamic_qtl_visualization_dir=$output_root"dynamic_qtl_hit_visualization/"


# Output directory containing gwas overlaps
gwas_overlap_dir=$output_root"gwas_overlap/"

# Output directory containing comparisons to time-step independent analysis
time_step_comparison_dir=$output_root"time_step_independent_comparison/"

# Output directory containing figures for manuscripts
manuscript_figures_dir=$output_root"manuscript_figures/"








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
            sh run_gaussian_dynamic_qtl.sh $dynamic_eqtl_input_file $output_file $model_version $permute $covariate_method $job_number $num_jobs
            fi
        done
    done
done






################### TO DO
# cell line overlap analysis visualization (previously 'visualize_cell_line_overlap_between_methods.R')
# Scatter comparing glm to glmm (previously in 'visualize_method_comparison.R')
# Visualize tissue-specific chromHMM enrichments (previously in 'make_plots_for_manuscript')
# Time step independent comparison visualizations (previously in 'visualize_time_step_independent_comparison.R')

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

model_version="glm"
covariate_method="pc1_5"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method

sh downstream_analysis_on_dynamic_eqtl_results.sh $model_version $covariate_method $num_jobs $parameter_string $dynamic_eqtl_input_file $qtl_results_dir $qtl_pvalue_distribution_visualization_dir $cell_line_overlap_analysis_dir $genotype_file $time_step_independent_stem $chrom_hmm_input_dir $tissue_specific_chrom_hmm_enrichment_dir $time_step_independent_comparison_dir $gsea_data_dir $gencode_file $gene_set_enrichment_dir















if false; then
sh multiple_testing_correction_and_visualization.sh $real_file $null_file $parameter_string $qtl_results_dir $num_jobs $dynamic_qtl_visualization_dir $input_data_file $model_version $covariate_method
fi



##########################################
# Step 4: Temporal QTL enrichment analysis
##########################################

covariate_method="pc1_5"

# All variant-gene pairs tested and their corresponding pvalues (using GLM)
gaussian_glm_qtl_results_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glm_covariate_method_"$covariate_method"_permute_False_results.txt"
# All significant variant gene pairs (eFDR <= .01) and their corresponding pvalues (using GLM)
gaussian_glm_qtl_significant_results=$qtl_results_dir$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glm_covariate_method_"$covariate_method"_efdr_.01_significant.txt"
# All significant genes (eFDR <= .01) and their top variant and their corresponding pvalues (using GLM)
gaussian_glm_qtl_significant_gene_results=$qtl_results_dir$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glm_covariate_method_"$covariate_method"_efdr_.01_significant_egenes.txt"
# All variant-gene pairs tested and their corresponding pvalues (using GLMM)
gaussian_glmm_qtl_results_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glmm_covariate_method_"$covariate_method"_permute_False_results.txt"
# All significant variant gene pairs (eFDR <= .01) and their corresponding pvalues (using GLMM)
gaussian_glmm_qtl_significant_results=$qtl_results_dir$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glmm_covariate_method_"$covariate_method"_efdr_.01_significant.txt"
# All significant genes (eFDR <= .01) and their top variant and their corresponding pvalues (using GLMM)
gaussian_glmm_qtl_significant_gene_results=$qtl_results_dir$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glmm_covariate_method_"$covariate_method"_efdr_.01_significant_egenes.txt"


# Dynamic qtl results file from NB method
if [ $covariate_method == "none" ]; then
    nb_dynamic_qtl_file=$nb_dynamic_qtl_dir"te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_dosage_covariate_method_none_permutation_scheme_none_permute_False_merged_dynamic_qtl_results.txt"
else
    nb_dynamic_qtl_file=$nb_dynamic_qtl_dir"te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_dosage_covariate_method_cell_line_fixed_effect_"$covariate_method"Xtime_permutation_scheme_none_permute_False_merged_dynamic_qtl_results.txt"
fi


output_stem=$qtl_visualization_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_covariate_method_"$covariate_method
if false; then
sbatch temporal_qtl_enrichment_analysis.sh $gaussian_glm_qtl_results_file $gaussian_glm_qtl_significant_results $gaussian_glm_qtl_significant_gene_results $gaussian_glmm_qtl_results_file $gaussian_glmm_qtl_significant_results $gaussian_glmm_qtl_significant_gene_results $nb_dynamic_qtl_file $genotype_file $output_stem $t_0_eqtl_results_file
fi

if false; then
Rscript visualize_cell_line_overlap_between_methods.R $qtl_visualization_dir "200"
fi





##########################################
# Step 5: ChromHMM enrichment analysis
##########################################

num_permutations="1000"
model_version="glmm_time"
covariate_method="pc1_6"


if false; then
top_genes_arr=( "50" "75" "100" "125" "150" "175" "200" "300" "400" "500" "1000" "1500" "2000" "2500" "3000" "4000" "5000" "6000" )
for top_genes in "${top_genes_arr[@]}"
do
	parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
	real_dynamic_qtl_results_file=$qtl_results_dir$parameter_string"_permute_False_results.txt"
	echo $parameter_string
	sbatch chrom_hmm_enrichment_analysis_cross_top_genes.sh $parameter_string$top_genes"_" $num_permutations $chrom_hmm_input_dir $real_dynamic_qtl_results_file $time_step_independent_stem $chrom_hmm_enrichment_directory $top_genes
done
fi


base_parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version
if false; then
Rscript visualize_chrom_hmm_enrichment_cross_top_genes.R $base_parameter_string $num_permutations $chrom_hmm_enrichment_directory
fi






####################################################
# Step 6: Tissue specific chrom-hmm enrichments
####################################################
model_version="glm"
covariate_method="pc1_5"
num_permutations="1000"
num_genes="550"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method


variant_gene_pairs_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_False_results.txt"


cluster_assignment_stem=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_sorted_variant_gene_pairs_"
##########################################################################
##########################################################################
cluster_assignment="dynamic_qtl_predicted_means_0.5"
##########################################################################
##########################################################################
if false; then
sh tissue_specific_chrom_hmm_enrichment_analysis.sh $parameter_string"_"$num_genes"_genes_" $num_permutations $chrom_hmm_input_dir $variant_gene_pairs_file $time_step_independent_stem $tissue_specific_chrom_hmm_enrichment_dir $cluster_assignment_stem $cluster_assignment $covariate_method $num_genes
fi


####################################################
# Step 7: Gene set enrichment analysis
####################################################
model_version="glm"
covariate_method="none"
num_genes="550"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method

variant_gene_pairs_file=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_genes_organized_input_data.txt"
cluster_assignment_stem=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_sorted_variant_gene_pairs_"
cluster_assignment="manual_annotations"
if false; then
sbatch gene_set_enrichment_analysis.sh $parameter_string $variant_gene_pairs_file $gencode_file $gene_set_enrichment_dir $gsea_data_dir $cluster_assignment_stem $cluster_assignment $time_step_independent_stem
fi



####################################################
# Step 8: Gene set enrichment analysis of cardiomyopathy genes
####################################################
model_version="glm"
covariate_method="none"
num_genes="550"

parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
variant_gene_pairs_file=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_genes_organized_input_data.txt"
all_hits_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_False_results.txt"
if false; then
python cardiomyopathy_gene_set_enrichment_analysis.py $cardiomyopathy_gene_list $variant_gene_pairs_file $all_hits_file $gencode_file $cardiomyopathy_gene_set_enrichment_dir $parameter_string
fi



####################################################
# Step 9:GWAS Overlap
####################################################
model_version="glm_quadratic"
covariate_method="pc1_5"
num_genes="693"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
variant_gene_pairs_file=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_variant_gene_pairs_organized_input_data.txt"
all_hits_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_False_results.txt"
if false; then
python compute_gwas_overlaps.py $gwas_catalog $gwas_overlap_dir $variant_gene_pairs_file $parameter_string $all_hits_file
fi


model_version="glm_quadratic"
covariate_method="pc1_5"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
num_genes="693"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
variant_gene_pairs_file=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_variant_gene_pairs_organized_input_data.txt"
all_hits_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_False_results.txt"
threshold="5e-8"
if false; then
python compute_gtex_gwas_overlaps.py $gtex_gwas_hits_dir $gwas_overlap_dir $variant_gene_pairs_file $parameter_string $all_hits_file $threshold
fi


####################################################
# Step 9: Time Step independent comparison file
####################################################
model_version="glm_quadratic"
covariate_method="pc1_5"
num_genes="693"
cluster_assignment="dynamic_qtl_quadratic_predicted_means_1.0"
# "dynamic_qtl_quadratic_predicted_means_1.0"
# "dynamic_qtl_predicted_means_1.0"


parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method


all_hits_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_False_results.txt"
if false; then
sh time_step_independent_comparison.sh $parameter_string"_top_"$num_genes"_" $time_step_independent_stem $time_step_comparison_dir $cluster_assignment $all_hits_file $num_genes
fi

if false; then
Rscript make_plots_for_manuscript.R $manuscript_figures_dir $dynamic_qtl_visualization_dir $time_step_comparison_dir $tissue_specific_chrom_hmm_enrichment_dir
fi





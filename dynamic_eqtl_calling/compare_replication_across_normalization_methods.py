import numpy as np 
import os
import sys
import pdb

def extract_egenes(egene_file):
	f = open(egene_file)
	egenes = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		rs_id = data[0]
		gene_id = data[1]
		egenes[rs_id + '_' + gene_id] = 0
	f.close()
	return egenes

def compute_replication_rate(egenes, rep1_file, pvalue_threshold):
	f = open(rep1_file)
	overlaps = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		rs_id = data[0]
		gene_id = data[1]
		pvalue = float(data[-1])
		test_name = rs_id + '_' + gene_id
		if test_name not in egenes:
			continue
		if pvalue < pvalue_threshold:
			overlaps[test_name] = 1
	f.close()
	fraction = float(len(overlaps))/len(egenes)
	return fraction


def check_replication(egene_file, rep1_file, rep2_file, nominal_pvalue_thresholds, output_file, transformation_type_1, transformation_type_2):
	egenes = extract_egenes(egene_file)
	t = open(output_file, 'w')
	t.write('pvalue_threshold\ttransformation_type\treplication_rate\n')

	for pvalue_threshold in nominal_pvalue_thresholds:
		rep_rate1 = compute_replication_rate(egenes, rep1_file, pvalue_threshold)
		t.write(str(pvalue_threshold) + '\t' + transformation_type_1 + '\t' + str(rep_rate1) + '\n')	
		rep_rate2 = compute_replication_rate(egenes, rep2_file, pvalue_threshold)
		t.write(str(pvalue_threshold) + '\t' + transformation_type_2 + '\t' + str(rep_rate2) + '\n')	
	t.close()




output_dir = sys.argv[1]




# Root directory for each normalization method:
rpkm_root="/project2/gilad/bstrober/ipsc_differentiation_19_lines/gaussian_dynamic_qtl_pipelines_v2/qtl_results/"
log_rpkm_root="/project2/gilad/bstrober/ipsc_differentiation_19_lines/gaussian_log_dynamic_eqtl_pipelines/qtl_results/"
inverse_normal_root="/project2/gilad/bstrober/ipsc_differentiation_19_lines/gaussian_gaussian_projected_dynamic_eqtl_pipelines/qtl_results/"


# File names
linear_dynamic_eqtl_egene_file_name = 'gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_efdr_.05_significant_egenes.txt'
nonlinear_dynamic_eqtl_egene_file_name = 'gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_efdr_.05_significant_egenes.txt'
linear_dynamic_eqtl_all_tests_file_name = 'gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_permute_False_results.txt'
nonlinear_dynamic_eqtl_all_tests_file_name = 'gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_permute_False_results.txt'


nominal_pvalue_thresholds = [.1, .01, .001, .0001]



# Version 1: Linear dynamic eQTLs 
# see if rpkm egenes replicate nominally in other two data sets
egene_file = rpkm_root + linear_dynamic_eqtl_egene_file_name
rep1_file = log_rpkm_root + linear_dynamic_eqtl_all_tests_file_name
rep2_file = inverse_normal_root + linear_dynamic_eqtl_all_tests_file_name

output_file = output_dir + 'replicatation_of_rpkm_linear_dynamic_eqtls_in_other_data_sets.txt'
check_replication(egene_file, rep1_file, rep2_file, nominal_pvalue_thresholds, output_file, 'log', 'inverse_normal')



# Version 1: Linear dynamic eQTLs 
# see if log(rpkm) egenes replicate nominally in other two data sets
egene_file = log_rpkm_root + linear_dynamic_eqtl_egene_file_name
rep1_file = rpkm_root + linear_dynamic_eqtl_all_tests_file_name
rep2_file = inverse_normal_root + linear_dynamic_eqtl_all_tests_file_name

output_file = output_dir + 'replicatation_of_log_rpkm_linear_dynamic_eqtls_in_other_data_sets.txt'
check_replication(egene_file, rep1_file, rep2_file, nominal_pvalue_thresholds, output_file, 'rpkm', 'inverse_normal')


# Version 1: Linear dynamic eQTLs 
# see if inverse_normal egenes replicate nominally in other two data sets
egene_file = inverse_normal_root + linear_dynamic_eqtl_egene_file_name
rep1_file = rpkm_root + linear_dynamic_eqtl_all_tests_file_name
rep2_file = log_rpkm_root + linear_dynamic_eqtl_all_tests_file_name

output_file = output_dir + 'replicatation_of_inverse_normal_linear_dynamic_eqtls_in_other_data_sets.txt'
check_replication(egene_file, rep1_file, rep2_file, nominal_pvalue_thresholds, output_file, 'rpkm', 'log')

import numpy as np 
import os
import sys
import pdb


def glm_extraction(gaussian_glm_qtl_results_file):
    f = open(gaussian_glm_qtl_results_file)
    pairs = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        rs_id = data[0]
        ensamble_id = data[1]
        pvalue = data[-1]
        pairs[rs_id + '_' + ensamble_id] = pvalue
    return pairs


def nb_extraction(nb_dynamic_qtl_file):
    f = open(nb_dynamic_qtl_file)
    head_count = 0
    pairs = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        rs_id = data[2]
        ensamble_id = data[5]
        pvalue = data[-3]
        pairs[rs_id + '_' + ensamble_id] = pvalue
    return pairs


qtl_results_dir = sys.argv[1]
visualization_input_dir = sys.argv[2]

output_file = visualization_input_dir + 'gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_covariate_method_pc1_5_merge_glm_glmm.txt'

gaussian_glm_qtl_results_file = qtl_results_dir + 'gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_permute_False_results.txt'
gaussian_glmm_qtl_results_file = qtl_results_dir + 'gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glmm_covariate_method_pc1_5_permute_False_results.txt'

glm_scores = glm_extraction(gaussian_glm_qtl_results_file)

glmm_scores = glm_extraction(gaussian_glmm_qtl_results_file)


t = open(output_file,'w')
t.write('glm_pvalue\tglmm_pvalue\n')

for key in glm_scores.keys():
    t.write(glm_scores[key] + '\t' + glmm_scores[key] + '\n')
t.close()
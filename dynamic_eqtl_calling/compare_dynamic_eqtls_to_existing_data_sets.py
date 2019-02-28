import numpy as np 
import os
import sys
import pdb




def extract_number_of_genes(sig_egene_file):
    f = open(sig_egene_file)
    counter = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        counter = counter + 1
    f.close()
    return counter



def extract_variant_gene_pair_time_step_info_using_dynamic_qtl_predicted_means(variant_gene_pairs_file, cluster_assignment_thresh):
    dicti = {}
    f = open(variant_gene_pairs_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        rs_id = data[0]
        ensamble_id = data[1]
        coefs = np.asarray(data[2].split(',')).astype(float)
        genotype_coef = coefs[1]
        interaction_coef = coefs[-1]
        predicted_time_0_genotype_0 = genotype_coef*0.0
        predicted_time_0_genotype_2 = genotype_coef*2.0
        predicted_time_15_genotype_0 = 15.0*0.0*interaction_coef + genotype_coef*0.0
        predicted_time_15_genotype_2 = 15.0*2.0*interaction_coef + genotype_coef*2.0

        t0 = predicted_time_0_genotype_0 - predicted_time_0_genotype_2
        t15 = predicted_time_15_genotype_0 - predicted_time_15_genotype_2
        test_name = rs_id + '_' + ensamble_id

        if t0 > 0 and t15 > 0:
            if abs(t0) > abs(t15):
                class_name = 'early'
            else:
                class_name = 'late'
        elif t0 <= 0 and t15 <= 0:
            if abs(t0) > abs(t15):
                class_name = 'early'
            else:
                class_name = 'late'
        else:
            if abs(t0) >= cluster_assignment_thresh and abs(t15) >= cluster_assignment_thresh:
                class_name = 'change'
            elif abs(t0) >= abs(t15):
                class_name = 'early'
            else:
                class_name = 'late'
        dicti[test_name] = class_name
    return dicti

def extract_variant_gene_pair_time_step_info_using_dynamic_qtl_quadratic_predicted_means(variant_gene_pairs_file, cluster_assignment_thresh):
    dicti = {}
    f = open(variant_gene_pairs_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        rs_id = data[0]
        ensamble_id = data[1]
        coefs = np.asarray(data[2].split(',')).astype(float)
        genotype_coef = coefs[1]
        interaction_coef = coefs[-2]
        squared_interaction_coef = coefs[-1]

        predicted_time_0_genotype_0 = genotype_coef*0.0
        predicted_time_0_genotype_2 = genotype_coef*2.0
        predicted_time_15_genotype_0 = 15.0*0.0*interaction_coef + genotype_coef*0.0 + 15.0*15.0*0.0*squared_interaction_coef
        predicted_time_15_genotype_2 = 15.0*2.0*interaction_coef + genotype_coef*2.0 + 15.0*15.0*2.0*squared_interaction_coef

        predicted_time_7_half_genotype_0 = 7.5*0.0*interaction_coef + genotype_coef*0.0 + 7.5*7.5*0.0*squared_interaction_coef
        predicted_time_7_half_genotype_2 = 7.5*2.0*interaction_coef + genotype_coef*2.0 + 7.5*7.5*2.0*squared_interaction_coef

        t0 = predicted_time_0_genotype_0 - predicted_time_0_genotype_2
        t_half = predicted_time_7_half_genotype_0 - predicted_time_7_half_genotype_2
        t15 = predicted_time_15_genotype_0 - predicted_time_15_genotype_2
        test_name = rs_id + '_' + ensamble_id
    
        if abs(t_half) >= (abs(t15)) and abs(t_half) >= (abs(t0)):
            class_name = 'middle'
        elif t0 > 0 and t15 > 0:
            if abs(t0) > abs(t15):
                class_name = 'early'
            else:
                class_name = 'late'
        elif t0 <= 0 and t15 <= 0:
            if abs(t0) > abs(t15):
                class_name = 'early'
            else:
                class_name = 'late'
        else:
            if abs(t0) >= cluster_assignment_thresh and abs(t15) >= cluster_assignment_thresh:
                class_name = 'change'
            elif abs(t0) >= abs(t15):
                class_name = 'early'
            else:
                class_name = 'late'
        dicti[test_name] = class_name
    return dicti


def extract_ipsc_results(file_name):
	dicti = {}
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		data = line.split()
		ensamble_id = data[0].split('.')[0]
		rs_id = data[1].split('.')[0]
		pvalue = float(data[3])
		test_name = ensamble_id + '_' + rs_id
		if rs_id == '':
			continue
		if test_name in dicti:
			print('assumption error!!')
			pdb.set_trace()
		dicti[test_name] = pvalue
	f.close()
	return dicti

def extract_ipsc_cm_results(file_name):
	f = open(file_name)
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[0]
		rs_id = data[1]
		pvalue = float(data[2])
		test_name = ensamble_id + '_' + rs_id
		if test_name in dicti:
			print('assumption error')
			pdb.set_trace()
		dicti[test_name] = pvalue
	return dicti

def print_results(sig_egene_file, cluster_assignments, qtl_results, output_file):
	t = open(output_file, 'w')
	t.write('rs_id\tensamble_id\tdynamic_eqtl_pvalue\tdynamic_eqtl_class\teqtl_pvalue\n')
	f = open(sig_egene_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		rs_id = data[0]
		ensamble_id = data[1]
		dynamic_pvalue = data[3]
		test_name = ensamble_id + '_' + rs_id
		if test_name in qtl_results:
			t.write(rs_id + '\t' + ensamble_id + '\t' + dynamic_pvalue + '\t' + cluster_assignments[rs_id + '_' + ensamble_id] + '\t' + str(qtl_results[test_name]) + '\n')
	t.close()
	f.close()
##########################
# Command Line args
##########################
parameter_string = sys.argv[1]
sig_egene_file = sys.argv[2]
all_hits_file = sys.argv[3]
model_version = sys.argv[4]
threshold = float(sys.argv[5])
output_dir = sys.argv[6]
ipsc_eqtl_file = sys.argv[7]
cm_eqtl_file = sys.argv[8]

if model_version == 'glm' or model_version == 'glmm':
	cluster_assignments = extract_variant_gene_pair_time_step_info_using_dynamic_qtl_predicted_means(all_hits_file, threshold)
elif model_version == 'glm_quadratic':
	cluster_assignments = extract_variant_gene_pair_time_step_info_using_dynamic_qtl_quadratic_predicted_means(all_hits_file, threshold)

ipsc_results = extract_ipsc_results(ipsc_eqtl_file)
ipsc_cm_results = extract_ipsc_cm_results(cm_eqtl_file)


ipsc_output_file = output_dir + parameter_string + '_' + str(threshold) + '_compare_dynamic_qtls_to_ipsc.txt'
print_results(sig_egene_file, cluster_assignments, ipsc_results, ipsc_output_file)

ipsc_cm_output_file = output_dir + parameter_string + '_' + str(threshold) + '_compare_dynamic_qtls_to_ipsc_cm.txt'
print_results(sig_egene_file, cluster_assignments, ipsc_cm_results, ipsc_cm_output_file)


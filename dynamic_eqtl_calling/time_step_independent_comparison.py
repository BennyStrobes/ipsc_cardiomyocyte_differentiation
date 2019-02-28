import numpy as np
import os
import sys
import pdb
import math
import random

def get_variant_gene_pairs(all_hits_file, pvalue_mapping, num_genes):
    top_tests = get_top_genes(all_hits_file, num_genes)
    dicti = {}
    for test_name in top_tests.keys():
        pvalue = float(pvalue_mapping[test_name])
        dicti[test_name] = np.zeros(17)
        dicti[test_name][-1] = pvalue
    return dicti

def fill_in_dictionary_with_time_step_independent_file(variant_gene_pairs, time_step_independent_file, time_step):
    f = open(time_step_independent_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[1]
        rs_id = data[3]
        pvalue = float(data[-1])
        if rs_id + '_' + ensamble_id in variant_gene_pairs:
            variant_gene_pairs[rs_id + '_' + ensamble_id][time_step] = pvalue
    return variant_gene_pairs

def get_cluster_assignments(cluster_assignment_file):
    f = open(cluster_assignment_file)
    dicti = {}
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        dicti[data[1] + '_' + data[0]] = data[2]
    f.close()
    return dicti

def get_top_genes(file_name, num_genes):
    f = open(file_name)
    genes = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        rs_id = data[0]
        ensamble_id = data[1]
        pvalue = float(data[-1])
        if ensamble_id not in genes:
            genes[ensamble_id] = (rs_id, pvalue)
        else:
            old_tupler = genes[ensamble_id]
            old_pvalue = old_tupler[1]
            if pvalue < old_pvalue:
                genes[ensamble_id] = (rs_id, pvalue)
    f.close()
    arr = []
    for geney in genes.keys():
        rs_id = genes[geney][0]
        pvalue = genes[geney][1]
        test_name = rs_id + '_' + geney 
        arr.append((test_name, pvalue))
    arr.sort(key=lambda x: x[1])
    top_tests = {}
    for i, tupler in enumerate(arr):
        if i >= num_genes:
            continue
        top_tests[tupler[0]] = 1
    print(len(top_tests))
    return top_tests

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

def extract_number_of_genes(sig_egene_file):
    f = open(sig_egene_file)
    counter = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        counter = counter + 1
    f.close()
    return counter

def make_all_tests_comparison_file(all_hits_file, sig_egene_file, time_step_independent_stem, dynamic_egenes_comparison_file, model_version, threshold, pvalue_mapping):
    # Extract number of significant genes
    num_genes = extract_number_of_genes(sig_egene_file)
    variant_gene_pairs = get_variant_gene_pairs(all_hits_file, pvalue_mapping, num_genes)
    if model_version == 'glm' or model_version == 'glmm':
        cluster_assignments = extract_variant_gene_pair_time_step_info_using_dynamic_qtl_predicted_means(all_hits_file, threshold)
    elif model_version == 'glm_quadratic':
        cluster_assignments = extract_variant_gene_pair_time_step_info_using_dynamic_qtl_quadratic_predicted_means(all_hits_file, threshold)

    names = []
    names.append('test_name')
    for time_step in range(16):
        print(time_step)
        time_step_independent_file =  time_step_independent_stem + str(time_step) + '_eqtl_results.txt'
        variant_gene_pairs = fill_in_dictionary_with_time_step_independent_file(variant_gene_pairs, time_step_independent_file, time_step)
        names.append('time_step_' + str(time_step))
    names.append('dynamic')
    t = open(dynamic_egenes_comparison_file, 'w')
    t.write('\t'.join(names) + '\t' + 'cluster_assignment' + '\n')
    for test_name in variant_gene_pairs.keys():
        t.write(test_name + '\t')
        pvalues = variant_gene_pairs[test_name].astype(str)
        t.write('\t'.join(pvalues) + '\t' + cluster_assignments[test_name] + '\n')
    t.close()


# Return the bin number corresponding to this distance
def get_distance_bin(distance, distance_bin_size):
    return int(math.floor(distance/distance_bin_size))


# Return the bin number corresponding to this distance
def get_maf_bin(maf, maf_bin_size):
    return int(math.floor(maf/maf_bin_size))

def get_background_variant_gene_pairs(input_dicti, input_file):
    distance_bin_size = 10000
    maf_bin_size = .05
    eqtl_distance = 50000
    ####################
    # Initialize object
    ####################
    background_qtls = []
    # number of bins needed for maf and distance
    num_distance_bins = int(math.ceil(eqtl_distance/distance_bin_size + 1))
    num_maf_bins = int(math.ceil(.5/maf_bin_size + 1))
    maf_mapping = {}
    distance_mapping = {}
    # Add each possible bin
    for distance_bin in range(num_distance_bins):
        background_qtls.append([])
        for maf_bin in range(num_maf_bins):
            background_qtls[distance_bin].append([])
    output_dicti = {}
    f = open(input_file)
    head_count = 0
    background = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[1]
        rs_id = data[3]
        test_name = rs_id + '_' + ensamble_id
        distance = abs(float(data[2]) - float(data[4]))
        maf = float(data[5])
        if test_name in input_dicti:
            maf_mapping[test_name] = maf
            distance_mapping[test_name] = distance
        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(distance, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)
        background_qtls[distance_bin][maf_bin].append(test_name)
    f.close()

    output_dicti = {}
    for test_name in input_dicti.keys():
        arr = input_dicti[test_name]
        distance = distance_mapping[test_name]
        maf = maf_mapping[test_name]
        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(distance, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)

        tests_to_choose_from = background_qtls[distance_bin][maf_bin]
        if len(tests_to_choose_from) < 4:
            print(len(tests_to_choose_from))
            print('short supply')
            pdb.set_trace()
        working = False
        while working == False:
            new_test = random.choice(tests_to_choose_from)
            if new_test not in output_dicti:
                output_dicti[new_test] = arr
                working = True
    return output_dicti
     

def make_background_tests_comparison_file(all_hits_file, sig_egene_file, time_step_independent_stem, all_tests_comparison_file, pvalue_mapping):
    num_genes = extract_number_of_genes(sig_egene_file)
    variant_gene_pairs_real = get_variant_gene_pairs(all_hits_file, pvalue_mapping, num_genes)
    variant_gene_pairs = get_background_variant_gene_pairs(variant_gene_pairs_real, time_step_independent_stem + '0_eqtl_results.txt' )
    names = []
    names.append('test_name')
    for time_step in range(16):
        print(time_step)
        time_step_independent_file = time_step_independent_stem + str(time_step) + '_eqtl_results.txt'
        variant_gene_pairs = fill_in_dictionary_with_time_step_independent_file(variant_gene_pairs, time_step_independent_file, time_step)
        names.append('time_step_' + str(time_step))
    names.append('dynamic')
    t = open(all_tests_comparison_file, 'w')
    t.write('\t'.join(names) + '\n')
    for test_name in variant_gene_pairs.keys():
        t.write(test_name + '\t')
        pvalues = variant_gene_pairs[test_name].astype(str)
        t.write('\t'.join(pvalues) + '\n')
    t.close()


def create_pvalue_mapping(all_hits_file):
    dicti = {}
    f = open(all_hits_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        rs_id = data[0]
        ensamble_id = data[1]
        pvalue = data[3]
        dicti[rs_id + '_' + ensamble_id] = pvalue
    return dicti


# Output files
dynamic_egenes_comparison_file = sys.argv[1]
dynamic_egenes_background_comparison_file = sys.argv[2]
# Input files
time_step_independent_stem = sys.argv[3]
threshold = float(sys.argv[4])
model_version = sys.argv[5]
all_hits_file = sys.argv[6]
sig_egene_file = sys.argv[7]


pvalue_mapping = create_pvalue_mapping(all_hits_file)


make_all_tests_comparison_file(all_hits_file, sig_egene_file, time_step_independent_stem, dynamic_egenes_comparison_file, model_version, threshold, pvalue_mapping)

make_background_tests_comparison_file(all_hits_file, sig_egene_file, time_step_independent_stem, dynamic_egenes_background_comparison_file, pvalue_mapping)
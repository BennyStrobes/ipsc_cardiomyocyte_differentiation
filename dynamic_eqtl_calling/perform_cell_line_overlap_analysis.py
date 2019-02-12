import numpy as np 
import os
import sys
import pdb
import math
import random


# Create dictionary list of all tested variants
def extract_variants_from_dynamic_qtl_file(dynamic_qtl_all_hits_file):
    f = open(dynamic_qtl_all_hits_file)
    head_count = 0
    variants = {} 
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        rs_id = data[0]
        variants[rs_id] = 0.0
    return variants

# Create dictionary list of all tested variants
def extract_variants_from_standard_qtl_file(dynamic_qtl_all_hits_file):
    f = open(dynamic_qtl_all_hits_file)
    head_count = 0
    variants = {} 
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        rs_id = data[3]
        variants[rs_id] = 0.0
    return variants


def extract_top_n_variants_from_dynamic_qtl_file(dynamic_qtl_all_hits_file, nn):
    f = open(dynamic_qtl_all_hits_file)
    head_count = 0
    variants = {}
    genes = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        rs_id = data[0]
        pvalue = float(data[-1])
        ensamble_id = data[1]
        if ensamble_id not in genes:
            genes[ensamble_id] = (pvalue, rs_id)
        else: # seen gene before
            old_tupler = genes[ensamble_id]
            old_pvalue = old_tupler[0]
            old_rs_id = old_tupler[1]
            if old_pvalue < pvalue:
                genes[ensamble_id] = old_tupler
            else:
                genes[ensamble_id] = (pvalue, rs_id)
    listy = []
    for ensamble_id in genes.keys():
        listy.append(genes[ensamble_id])
    sorted_listy = sorted(listy, key=lambda x: x[0])
    for i, tupler in enumerate(sorted_listy):
        variants[tupler[1]] = 0
        if len(variants) == nn:
            return variants

def extract_top_n_variants_from_standard_qtl_file(dynamic_qtl_all_hits_file, nn):
    f = open(dynamic_qtl_all_hits_file)
    head_count = 0
    variants = {}
    genes = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        rs_id = data[3]
        pvalue = float(data[-1])
        ensamble_id = data[1]
        if ensamble_id not in genes:
            genes[ensamble_id] = (pvalue, rs_id)
        else: # seen gene before
            old_tupler = genes[ensamble_id]
            old_pvalue = old_tupler[0]
            old_rs_id = old_tupler[1]
            if old_pvalue < pvalue:
                genes[ensamble_id] = old_tupler
            else:
                genes[ensamble_id] = (pvalue, rs_id)
    listy = []
    for ensamble_id in genes.keys():
        listy.append(genes[ensamble_id])
    sorted_listy = sorted(listy, key=lambda x: x[0])
    for i, tupler in enumerate(sorted_listy):
        variants[tupler[1]] = 0
        if len(variants) == nn:
            return variants


# Create dictionary list of all tested variants
def extract_sig_variants_from_dynamic_qtl_file(file_name):
    f = open(file_name)
    head_count = 0
    variants = {} 
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        rs_id = data[2]
        variants[rs_id] = 0.0
    return variants

# Compute minor allele frequency
def get_maf(dosages):
    af = np.sum(dosages)/(2*len(dosages))
    if af > .5:
        maf = 1.0 - af
    else:
        maf = af

    if maf > .5 or maf < .1:
        print(maf)
    return maf

# Create mapping from variants to MAF (only have to create for tested_variants cause includes sig_variants)
def create_mapping_from_variant_to_maf(tested_variants, genotype_file):
    # Initialize mapping
    variant_to_maf = {}
    # Stream genotype file
    f = open(genotype_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#CHROM'):
            cell_lines = np.asarray(data[9:])
        if line.startswith('#'):
            continue
        rs_id = data[2]
        if rs_id not in tested_variants:  # Did not test this variant so don't care about it
            continue
        dosages = np.asarray(data[9:]).astype(float)
        maf = get_maf(dosages)
        genotype = np.round(dosages)
        variant_to_maf[rs_id] = (maf, genotype, dosages)
    return variant_to_maf, cell_lines


def compute_overlap_matrix(sig_variants, variant_to_maf, cell_line_names):
    num_lines = len(cell_line_names)
    overlap_matrix = np.zeros((num_lines, num_lines))
    for rs_id in sig_variants.keys():
        genotype = variant_to_maf[rs_id][1]
        for ii in range(num_lines):
            for jj in range(num_lines):
                if ii == jj:
                    continue
                if genotype[ii] == genotype[jj]:
                    overlap_matrix[ii, jj] = overlap_matrix[ii,jj] + 1
    return overlap_matrix


# Return the bin number corresponding to this distance
def get_maf_bin(maf, maf_bin_size):
    return int(math.floor(maf/maf_bin_size))

def get_sample_group_size(genotype):
    aa = len(np.where(genotype==0.0)[0])
    bb = len(np.where(genotype==1.0)[0])
    cc = len(np.where(genotype==2.0)[0])
    return max(aa,bb,cc)

def get_middle_sample_group_size(genotype):
    aa = len(np.where(genotype==0.0)[0])
    bb = len(np.where(genotype==1.0)[0])
    cc = len(np.where(genotype==2.0)[0])
    return sorted([aa,bb,cc])[1]

def make_background_distribution(variant_to_maf, maf_bin_size, sample_group_bin_size):
    ####################
    # Initialize object
    ####################
    background_qtls = []
    # number of bins needed for maf and distance
    num_maf_bins = int(math.ceil(.5/maf_bin_size + 1))

    num_group_bins = int(math.ceil(19/sample_group_bin_size + 1))
    # Add each possible bin
    for maf_bin in range(num_maf_bins):
        background_qtls.append([])
        for group_bin in range(num_group_bins):
            background_qtls[maf_bin].append([])
            for group_bin2 in range(num_group_bins):
                background_qtls[maf_bin][group_bin].append([])
    for rs_id in variant_to_maf.keys():
        maf = variant_to_maf[rs_id][0]
        genotype = variant_to_maf[rs_id][1]
        sample_group_size = get_sample_group_size(genotype)
        sample_group_size2 = get_middle_sample_group_size(genotype)
        maf_bin = get_maf_bin(maf, maf_bin_size)
        background_qtls[maf_bin][sample_group_size][sample_group_size2].append(rs_id)
    return background_qtls

def get_background_variants(background_distribution, sig_variants, variant_to_maf, maf_bin_size):
    background_variants_perm = {}
    for variant in sig_variants.keys():
        maf = variant_to_maf[variant][0]
        maf_bin = get_maf_bin(maf, maf_bin_size)
        genotype = variant_to_maf[variant][1]
        sample_group_size = get_sample_group_size(genotype)
        sample_group_size2 = get_middle_sample_group_size(genotype)
        available_variants = background_distribution[maf_bin][sample_group_size][sample_group_size2]
        converged = False
        while converged == False:
            new_variant = random.choice(available_variants)
            if new_variant not in background_variants_perm:
                background_variants_perm[new_variant] = 1
                converged = True
    return background_variants_perm

def compute_pvalues(overlap_matrix_real, overlap_matrices_perm):
    num_cell_lines = overlap_matrix_real.shape[0]
    pvalue_matrix = np.zeros((num_cell_lines, num_cell_lines))
    num_perms = len(overlap_matrices_perm)
    for ii in range(num_cell_lines):
        for jj in range(num_cell_lines):
            if ii == jj:
                continue
            counter = 0
            for itera in range(num_perms):
                if overlap_matrix_real[ii,jj] <  overlap_matrices_perm[itera][ii,jj]:
                    counter = counter + 1
            pvalue = float(counter)/num_perms
            pvalue_matrix[ii,jj] = pvalue
    return pvalue_matrix


def statistical_analysis(tested_variants, sig_variants, variant_to_maf, cell_line_names):
    # Produce overlap matrix for real data
    overlap_matrix_real = compute_overlap_matrix(sig_variants, variant_to_maf, cell_line_names)
    # Produce overlap matrix for permuted data
    num_perms = 10000
    maf_bin_size = .01
    sample_group_bin_size = 1
    overlap_matrices_perm = []
    background_distribution = make_background_distribution(variant_to_maf, maf_bin_size, sample_group_bin_size)
    for perm_num in range(num_perms):
        print(perm_num)
        # randomly select len(sig_variants) variants from tested_variants that are matched for MAF
        background_variants_perm = get_background_variants(background_distribution, sig_variants, variant_to_maf, maf_bin_size)
        overlap_matrix_perm = compute_overlap_matrix(background_variants_perm, variant_to_maf, cell_line_names)
        overlap_matrices_perm.append(overlap_matrix_perm)

    pvalue_matrix = compute_pvalues(overlap_matrix_real, overlap_matrices_perm)
    return pvalue_matrix

def compute_perm_overlap_matrix(tested_variants, sig_variants, variant_to_maf, cell_line_names, num_perms):
    maf_bin_size = .01
    sample_group_bin_size = 1
    overlap_matrices_perm = []
    background_distribution = make_background_distribution(variant_to_maf, maf_bin_size, sample_group_bin_size)
    overlap_matrix_perm = np.zeros((len(cell_line_names),len(cell_line_names)))
    for perm_num in range(num_perms):
        # randomly select len(sig_variants) variants from tested_variants that are matched for MAF
        background_variants_perm = get_background_variants(background_distribution, sig_variants, variant_to_maf, maf_bin_size)
        overlap_matrix_perm_curr = compute_overlap_matrix(background_variants_perm, variant_to_maf, cell_line_names)
        overlap_matrix_perm = overlap_matrix_perm + overlap_matrix_perm_curr
    return overlap_matrix_perm/(len(sig_variants)*num_perms)

def print_overlap_matrix(overlap_matrix, cell_line_names, output_file):
    t = open(output_file, 'w')
    t.write('cell_line\t' + '\t'.join(cell_line_names) + '\n')
    for i, cell_line_name in enumerate(cell_line_names):
        t.write(cell_line_name + '\t' + '\t'.join(overlap_matrix[i,:].astype(str)) + '\n')
    t.close()

# Create dictionary list of all time-step significant variants
def extract_sig_variants_from_time_step_independent_file(time_step_independent_stem):
    varz = {}
    for time_step in range(16):
        file_name = time_step_independent_stem + str(time_step) + '_efdr_thresh_.1_significant_egenes.txt'
        f = open(file_name)
        head_count = 0
        for line in f:
            line = line.rstrip()
            data = line.split()
            if head_count == 0:
                head_count = head_count + 1
                continue
            rs_id = data[3]
            varz[rs_id] = 0
        f.close()
    return varz

def print_dosage_matrix(sig_variants, variant_to_maf, output_file, cell_line_names):
    t = open(output_file, 'w')
    t.write('\t'.join(cell_line_names) + '\n')
    for variant in sig_variants.keys():
        dosage = variant_to_maf[variant][2].astype(str)
        t.write('\t'.join(dosage) + '\n')
    t.close()


qtl_all_hits_file = sys.argv[1]
genotype_file = sys.argv[2]
real_overlap_matrix_file = sys.argv[3]
perm_overlap_matrix_file = sys.argv[4]
num_hits = int(sys.argv[5])
version = sys.argv[6]


if version == 'dynamic_eqtl':
    # Create dictionary list of all tested variants
    tested_variants = extract_variants_from_dynamic_qtl_file(qtl_all_hits_file)
    # Extract top n variants
    top_nn_variants = extract_top_n_variants_from_dynamic_qtl_file(qtl_all_hits_file, num_hits)
elif version == 'standard_eqtl':
    # Create dictionary list of all tested variants
    tested_variants = extract_variants_from_standard_qtl_file(qtl_all_hits_file)
    # Extract top n variants
    top_nn_variants = extract_top_n_variants_from_standard_qtl_file(qtl_all_hits_file, num_hits)


# Create mapping from variants to:
### 1. MAF (only have to create for tested_variants cause includes sig_variants)
### 2. genotype vector
variant_to_maf, cell_line_names = create_mapping_from_variant_to_maf(tested_variants, genotype_file)


# For top num_Hits
# Compute matrix of dimension num_cell_lineXnum_cell_line where each element in the matrix represents the fraction of times those two cell lines were in the same genotype class
overlap_matrix_real_counts = compute_overlap_matrix(top_nn_variants, variant_to_maf, cell_line_names)
overlap_matrix_real = overlap_matrix_real_counts/len(top_nn_variants)

num_perms=1
overlap_matrix_perm = compute_perm_overlap_matrix(tested_variants, top_nn_variants, variant_to_maf, cell_line_names, num_perms)


print_overlap_matrix(overlap_matrix_real, cell_line_names, real_overlap_matrix_file)

print_overlap_matrix(overlap_matrix_perm, cell_line_names, perm_overlap_matrix_file)


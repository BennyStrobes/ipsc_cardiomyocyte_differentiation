import numpy as np
import os
import sys
import pdb

# First, get index of sample that has the largest qvalue less than qval_thresh
def get_sample_index(qvalue_file, qval_thresh):
    f = open(qvalue_file)
    min_diff = 1000000
    best_index = -5
    for counter, line in enumerate(f):
        line = line.rstrip()
        data = line.split()
        if counter == 0:
            continue
        qval = float(data[1])
        diff = abs(qval - qval_thresh)
        if qval < qval_thresh and diff < min_diff:
            min_diff = diff
            best_index = counter
    return best_index

        
def get_pvalz(file_name):
    aa = np.loadtxt(file_name,dtype=str,delimiter='\t')
    return np.asarray(aa[1:,-1]).astype(float)

# Get pvalue threshold corresponding to qval_thresh in null
def get_pval_thresh(efdr_file, fdr_thresh):
    f = open(efdr_file)
    head_count = 0
    prev_pval = -1.0
    f = open(efdr_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        nominal_pvalue = float(data[1])
        efdr = float(data[3])
        if efdr > fdr_thresh:
            return prev_pval
        else:
            prev_pval = nominal_pvalue
    print('fatal error with efdr: did not include enough top-snps param!')
    pdb.set_trace()

def filter_to_significant_results(real_file, significant_results_file, pval_thresh):
    f = open(real_file)
    t = open(significant_results_file, 'w')
    head_count = 0 
    for line in f:
        line = line.rstrip()
        data = line.split()
        pval = float(data[-1])
        if pval > pval_thresh:  # Not significant
            continue
        t.write(line + '\n')
    t.close()

def filter_to_egenes(significant_results_file, significant_gene_results_file):
    f = open(significant_results_file)
    t = open(significant_gene_results_file,'w')
    genes = {}
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        pval = float(data[-1])
        gene = data[1]
        if gene not in genes:
            genes[gene] = (pval, line)
        else:
            old_pval = genes[gene][0]
            old_line = genes[gene][1]
            if pval < old_pval:
                genes[gene] = (pval, line)
            else:
                genes[gene] = (old_pval,old_line)
    f.close()
    for gene in genes.keys():
        liner = genes[gene][1]
        t.write(liner + '\n')
    t.close()

efdr_file = sys.argv[1]
real_file = sys.argv[2]
significant_results_file = sys.argv[3]  #  Output file
significant_gene_results_file = sys.argv[4]  # Ouptut file
fdr_thresh = float(sys.argv[5])



# Get pvalue threshold corresponding to fdr_thresh in null
pval_thresh = get_pval_thresh(efdr_file, fdr_thresh)

# Create list of variant-gene pairs that have pvalue less than pvalue threshold
filter_to_significant_results(real_file, significant_results_file, pval_thresh)

# Create list of significant genes and their highest associated variant
filter_to_egenes(significant_results_file, significant_gene_results_file)
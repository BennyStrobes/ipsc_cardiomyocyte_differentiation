import numpy as np 
import os
import sys
import gzip
import pdb
from sklearn.decomposition import NMF
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("white")
sns.set(font_scale=1.4)


def visualize_heatmap(H, heatmap_file, column_labels,title):
    num_components = H.shape[0]
    row_labels = np.arange(1, num_components+1).astype(str)

    fig = plt.figure()
    sns.heatmap(H,cmap="Blues", cbar_kws={'label': 'Loadings'})
    plt.xlabel('Time Step (day)')
    plt.ylabel('Latent Factor')
    plt.title('Summary statistic: ' + title)
    plt.yticks(rotation=90)

    fig.savefig(heatmap_file)


def reorder_latent_factors(W,H):
    time_point_of_largest_loading_per_factor = np.argmax(H,axis=1)
    indexes = np.argsort(time_point_of_largest_loading_per_factor)
    return W[:,indexes], H[indexes,:]

def standardize_rows(input_mat):
    num_rows = input_mat.shape[0]
    for row_num in range(num_rows):
        standardized = (input_mat[row_num,:] - np.mean(input_mat[row_num,:]))/np.std(input_mat[row_num,:])
        #input_mat[row_num,:] = standardized - np.min(standardized)
        input_mat[row_num,:] = standardized
    input_mat = input_mat - np.min(input_mat)
    return input_mat

def intercept_shift_rows(input_mat):
    num_rows = input_mat.shape[0]
    for row_num in range(num_rows):
        standardized = (input_mat[row_num,:] - np.mean(input_mat[row_num,:]))
        input_mat[row_num,:] = standardized
    input_mat = input_mat - np.min(input_mat)
    return input_mat

def intercept_shift_columns(input_mat):
    num_colums = input_mat.shape[1]
    for col_num in range(num_colums):
        standardized = (input_mat[:,col_num] - np.mean(input_mat[:,col_num]))
        input_mat[:,col_num] = standardized
    input_mat = input_mat - np.min(input_mat)
    return input_mat

def run_matrix_factorization(input_mat, row_labels, column_labels, output_root, title, num_components, alpha_param):
    #input_mat = standardize_rows(input_mat)
    #input_mat = intercept_shift_rows(input_mat)
    #input_mat = intercept_shift_columns(input_mat)
    model = NMF(n_components=num_components, init='nndsvd', alpha=alpha_param, l1_ratio=1)
    W = model.fit_transform(input_mat)
    H = model.components_

    W, H = reorder_latent_factors(W,H)

    # Save loading matrix to output file
    loading_matrix_output_file = output_root + 'loading_matrix.txt'
    np.savetxt(loading_matrix_output_file, H, delimiter='\t',fmt="%s")

    # Make heatmap visualizing H
    heatmap_file = output_root + 'heatmap.png'
    visualize_heatmap(H, heatmap_file, column_labels,title)
    return W

# Extract list of all genes used as input to wasp
def extract_list_of_all_genes(cht_input_file):
    all_genes = {}
    f = open(cht_input_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        ensamble_id = data[6].split('_')[2]
        all_genes[ensamble_id] = 1
    return all_genes


# Extract list of all egenes found in WASP
def extract_list_of_egenes(variant_gene_pairs):
    egenes = {}
    for variant_gene in variant_gene_pairs:
        ensamble_id = variant_gene.split('_')[0]
        egenes[ensamble_id] = 1
    return egenes

# Extract array of length number of components. Where each element of array is a dictionary that contains genes in the test set for that componenent
def extract_test_set_genes(W, variant_gene_pairs, num_components, num_genes):
    test_sets = []
    num_rows = W.shape[0]
    for component_num in np.arange(num_components):
        genes = {}
        counter = 0
        loadings = W[:, component_num]
        top_n_indices = loadings.argsort()[-num_genes:][::-1]
        for index in top_n_indices:
            ensamble_id = variant_gene_pairs[index].split('_')[0]
            genes[ensamble_id] = 1
            counter = counter + 1
        test_sets.append(genes)
    return test_sets

def convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file):
    f = gzip.open(gencode_file)
    gene_symbol_hits = []
    for line in f:
        line = line.decode('utf-8').rstrip()
        data = line.split()
        if line.startswith('#'):
            continue
        line_ensamble_id = data[9].split('"')[1].split('.')[0]
        line_gene_symbol = data[17].split('"')[1]
        if line_ensamble_id in ensamble_hits:
            gene_symbol_hits.append(line_gene_symbol)
    return np.unique(gene_symbol_hits)

def sort_gsea(save_file, new_save_file):
    f = open(save_file)
    t = open(new_save_file,'w')
    pairs = []
    for i,line in enumerate(f):
        line = line.rstrip()
        data = line.split()
        if i < 4:
            continue
        pvalue = float(data[6])
        pairs.append((pvalue, line))
    sorted_pairs = sorted(pairs, key=lambda x: x[0])
    for pair in sorted_pairs:
        liner = pair[1]
        t.write(liner + '\n')
    t.close()

def sort_gsea_sig_only(save_file, new_save_file):
    f = open(save_file)
    t = open(new_save_file,'w')
    pairs = []
    for i,line in enumerate(f):
        line = line.rstrip()
        data = line.split()
        if i < 4:
            continue
        pvalue = float(data[6])
        pairs.append((pvalue, line))
    sorted_pairs = sorted(pairs, key=lambda x: x[0])
    for pair in sorted_pairs:
        liner = pair[1]
        if float(pair[1].split('\t')[7]) <= .05:
            t.write(liner + '\n')
    t.close()

def print_array(file_name, array):
    t = open(file_name,'w')
    for ele in array:
        t.write(ele + '\n')
    t.close()

def gene_set_enrichment(ensamble_hits, ensamble_background, output_root, cht_visualization_dir, gencode_file):
    gene_symbol_hits = convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file)
    gene_symbol_background = convert_from_ensamble_to_gene_symbol(ensamble_background, gencode_file)


    hits_file = output_root + 'temp_hit_genes.txt'
    background_file = output_root + 'temp_background_genes.txt'
    print_array(hits_file, gene_symbol_hits)
    print_array(background_file, gene_symbol_background)

    gene_sets = {'hallmark':'h.all.v5.1.symbols.gmt.txt'}

    for gene_set in gene_sets.keys():
        save_file = output_root + gene_set + '_gsea_output.txt'
        geneset_file = '/project2/gilad/bstrober/tools/tools/gsea/data/' + gene_sets[gene_set]
        os.system('gsea ' + hits_file + ' ' + background_file + ' ' + geneset_file + ' ' + save_file)
        new_save_file = output_root + gene_set + '_gsea_sorted_output.txt'
        sort_gsea(save_file, new_save_file)
        new_save_file_sig_only = output_root + gene_set + '_gsea_sorted_sig_only_output.txt'
        sort_gsea_sig_only(save_file, new_save_file_sig_only)
        os.system('rm ' + save_file)
    os.system('rm ' + hits_file)
    os.system('rm ' + background_file)



def gene_set_enrichment_shell(W, all_genes, egenes, variant_gene_pairs, output_root, cht_visualization_dir, gencode_file, num_genes):
    num_components = W.shape[1]
    # Extract array of length number of components. Where each element of array is a dictionary that contains genes in the test set for that componenent
    test_sets = extract_test_set_genes(W, variant_gene_pairs, num_components, num_genes)
    for component_num in np.arange(num_components):

        test_set = test_sets[component_num]
        print(len(test_set))
        # Run Enrichment for egenes only
        gene_set_enrichment(test_set, egenes, output_root + str(component_num) + '_egenes_', cht_visualization_dir, gencode_file)
        # Run Enrichment for all genes
        gene_set_enrichment(test_set, all_genes, output_root + str(component_num) + '_all_genes_', cht_visualization_dir, gencode_file)


#######################
# Command Line Args
#######################

parameter_string = sys.argv[1]
cht_output_dir = sys.argv[2]
matrix_factorization_dir = sys.argv[3]
target_regions_dir = sys.argv[4]
fdr = sys.argv[5]
pc_num = int(sys.argv[6])



# Extract list of all genes used as input to wasp
cht_input_file = target_regions_dir + 'target_regions_' + parameter_string + '_merged.txt'
all_genes = extract_list_of_all_genes(cht_input_file)


# Load in data
alpha_file = cht_output_dir + 'best_variant_per_egene_' + parameter_string + '_num_pc_' + str(pc_num) + '_fdr_' + fdr + '_alpha.txt'
beta_file = cht_output_dir + 'best_variant_per_egene_' + parameter_string + '_num_pc_' + str(pc_num) + '_fdr_' + fdr + '_beta.txt'
pvalue_file = cht_output_dir + 'best_variant_per_egene_' + parameter_string + '_num_pc_' + str(pc_num) +'_fdr_' + fdr + '_pvalues.txt'

alpha_full = np.loadtxt(alpha_file,dtype=str)
beta_full = np.loadtxt(beta_file, dtype=str)
pvalue_full = np.loadtxt(pvalue_file, dtype=str)

alpha = alpha_full[:, 1:].astype(float)
beta = beta_full[:, 1:].astype(float)
pvalue = pvalue_full[:, 1:].astype(float)

variant_gene_pairs = alpha_full[:,0]

allelic_fraction = np.abs((alpha/(alpha+beta)) - .5) 

alpha_params = [0, .25, .5, .75, 1, 1.5, 2, 5, 10]

log_pvalue = -np.log10(pvalue + 1e-17)

for num_components in np.arange(3, 6):
    for alpha_param in alpha_params:

        ## VERSION 1: allelic fraction
        version = 'allelic_fraction'
        output_root = matrix_factorization_dir + parameter_string + '_num_pc_' + str(pc_num) + '_fdr_' + fdr + '_' + version + '_factorization' + '_alpha_' + str(alpha_param) + '_' + str(num_components) + '_'
        W = run_matrix_factorization(allelic_fraction, variant_gene_pairs, np.arange(16).astype(str), output_root, version, num_components, alpha_param)


        ## VERSION 2: pvalue
        version = 'pvalue'
        output_root = matrix_factorization_dir + parameter_string + '_num_pc_' + str(pc_num) + '_fdr_' + fdr + '_' + version + '_factorization' + '_alpha_' + str(alpha_param) + '_' + str(num_components) + '_'
        W = run_matrix_factorization(pvalue, variant_gene_pairs, np.arange(16).astype(str), output_root, version, num_components, alpha_param)

        ## VERSION 3: -log10(pvalue)
        version = 'log_pvalue'
        output_root = matrix_factorization_dir + parameter_string + '_num_pc_' + str(pc_num) + '_fdr_' + fdr + '_' + version + '_factorization' + '_alpha_' + str(alpha_param) + '_' + str(num_components) + '_'
        W = run_matrix_factorization(log_pvalue, variant_gene_pairs, np.arange(16).astype(str), output_root, version, num_components, alpha_param)










#########################################
# Retired scripts (not currently used)
##########################################
# Extract list of all egenes found in WASP
#egenes = extract_list_of_egenes(variant_gene_pairs)



    #gene_set_enrichment_shell(W, all_genes, egenes, variant_gene_pairs, output_root, cht_visualization_dir, gencode_file, num_genes)

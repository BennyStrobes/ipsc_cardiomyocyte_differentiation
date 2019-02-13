import numpy as np 
import os
import sys
import pdb
import gzip
import scipy.stats


def extract_top_ensamble_ids(time_step_independent_file, significant_variant_gene_pairs_file, num_genes):
    ########################
    # Extract dictionary containing names of all tested genes
    all_genes = {}
    f = open(time_step_independent_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[1]
        all_genes[ensamble_id] = 1
    f.close()
    ########################
    # Extract dictionary list of genes that are signficant in this hits_version
    significant_genes = {}
    f = open(significant_variant_gene_pairs_file)
    head_count = 0
    used = {}
    listy = []
    for line in f:
        line = line.rstrip()
        data = line.split()
        ensamble_id = data[1]
        rs_id = data[0]
        pvalue = float(data[-1])
        if ensamble_id in used:
            pdb.set_trace()
        used[ensamble_id] = 1
        listy.append((ensamble_id, pvalue))
        #significant_genes[ensamble_id] = 1 
    f.close()
    listy.sort(key=lambda x: x[1])
    for i, tupler in enumerate(listy):
        if i < num_genes:
            significant_genes[tupler[0]] = 1
    return significant_genes, all_genes

def extract_ensamble_ids(time_step_independent_file, significant_variant_gene_pairs_file):
    ########################
    # Extract dictionary containing names of all tested genes
    all_genes = {}
    f = open(time_step_independent_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[1]
        all_genes[ensamble_id] = 1
    f.close()
    ########################
    # Extract dictionary list of genes that are signficant in this hits_version
    significant_genes = {}
    f = open(significant_variant_gene_pairs_file)
    head_count = 0
    used = {}
    listy = []
    for line in f:
        line = line.rstrip()
        data = line.split()
        ensamble_id = data[1]
        rs_id = data[0]
        pvalue = float(data[-1])
        significant_genes[ensamble_id] = 1
    f.close()
    return significant_genes, all_genes



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


def print_array(file_name, array):
    t = open(file_name,'w')
    for ele in array:
        t.write(ele + '\n')
    t.close()

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

def extract_variant_gene_pair_time_step_info_from_file(file_name):
    dicti = {}
    f = open(file_name)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        rs_id = data[0]
        gene_id = data[1]
        test_name = rs_id + '_' + gene_id
        cluster = data[2]
        dicti[test_name] = cluster 
    f.close()
    return dicti

#########################
# Command Line args
#########################
parameter_string = sys.argv[1]  # String that describes current version of the data
significant_variant_gene_pairs_file = sys.argv[2]  # File containing significant dynamic qtl hits
gencode_file = sys.argv[3]  # Gencode gene annotation file to be used to convert from ensamble id to gene-symbol id
gene_set_enrichment_directory = sys.argv[4]  # Output directory
gsea_data_dir = sys.argv[5]  # Directory containing all necessary gsea data
time_step_independent_stem = sys.argv[6]




# Results file for time step 0 (could have been any time step)
# Only use this to extract maf and distance to tss info
time_step_independent_file = time_step_independent_stem + '0_eqtl_results.txt'


ensamble_hits, ensamble_background = extract_ensamble_ids(time_step_independent_file, significant_variant_gene_pairs_file)


gene_symbol_hits = convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file)

gene_symbol_background = convert_from_ensamble_to_gene_symbol(ensamble_background, gencode_file)




hits_file = gene_set_enrichment_directory + parameter_string + '_hit_genes.txt'

background_file = gene_set_enrichment_directory + parameter_string + '_background_genes.txt'


print_array(hits_file, gene_symbol_hits)
print_array(background_file, gene_symbol_background)
#np.savetxt(hits_file, gene_symbol_hits,fmt="%s",delimiter="\n")
#np.savetxt(background_file, gene_symbol_background,fmt="%s",delimiter="\n")

#genesets = ['h.all.v5.1.symbols.gmt.txt', 'c2.cp.biocarta.v5.1.symbols.gmt.txt', 'c2.cp.kegg.v5.1.symbols.gmt.txt', 'c2.cp.reactome.v5.1.symbols.gmt.txt']

#names = ['hallmark', 'biocarta', 'kegg', 'reactome']

genesets = ['h.all.v5.1.symbols.gmt.txt']

names = ['hallmark']


for i, val in enumerate(genesets):
    name = names[i]
    geneset_file = gsea_data_dir + val
    save_file = gene_set_enrichment_directory + parameter_string + '_' + name + '_gsea_output.txt'
    #geneset_file = '/project2/gilad/bstrober/tools/tools/gsea/data/' + 'c2.cp.biocarta.v5.1.symbols.gmt.txt'
    os.system('gsea ' + hits_file + ' ' + background_file + ' ' + geneset_file + ' ' + save_file)


    new_save_file = gene_set_enrichment_directory + parameter_string + '_' + name + '_gsea_sorted_output.txt'
    sort_gsea(save_file, new_save_file)
    # Remove un-sorted file
    os.system('rm ' + save_file)

# Remove some unnecessary files
os.system('rm ' + hits_file)
os.system('rm ' + background_file)
























num_genes=200
parameter_string = parameter_string + '_top_' + str(num_genes) + '_genes'

# Results file for time step 0 (could have been any time step)
# Only use this to extract maf and distance to tss info
time_step_independent_file = time_step_independent_stem + '0_eqtl_results.txt'


ensamble_hits, ensamble_background = extract_top_ensamble_ids(time_step_independent_file, significant_variant_gene_pairs_file, num_genes)


gene_symbol_hits = convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file)

gene_symbol_background = convert_from_ensamble_to_gene_symbol(ensamble_background, gencode_file)




hits_file = gene_set_enrichment_directory + parameter_string + '_hit_genes.txt'

background_file = gene_set_enrichment_directory + parameter_string + '_background_genes.txt'


print_array(hits_file, gene_symbol_hits)
print_array(background_file, gene_symbol_background)
#np.savetxt(hits_file, gene_symbol_hits,fmt="%s",delimiter="\n")
#np.savetxt(background_file, gene_symbol_background,fmt="%s",delimiter="\n")

#genesets = ['h.all.v5.1.symbols.gmt.txt', 'c2.cp.biocarta.v5.1.symbols.gmt.txt', 'c2.cp.kegg.v5.1.symbols.gmt.txt', 'c2.cp.reactome.v5.1.symbols.gmt.txt']

#names = ['hallmark', 'biocarta', 'kegg', 'reactome']

genesets = ['h.all.v5.1.symbols.gmt.txt']

names = ['hallmark']


for i, val in enumerate(genesets):
    name = names[i]
    geneset_file = gsea_data_dir + val
    save_file = gene_set_enrichment_directory + parameter_string + '_' + name + '_gsea_output.txt'
    #geneset_file = '/project2/gilad/bstrober/tools/tools/gsea/data/' + 'c2.cp.biocarta.v5.1.symbols.gmt.txt'
    os.system('gsea ' + hits_file + ' ' + background_file + ' ' + geneset_file + ' ' + save_file)


    new_save_file = gene_set_enrichment_directory + parameter_string + '_' + name + '_gsea_sorted_output.txt'
    sort_gsea(save_file, new_save_file)
    # Remove un-sorted file
    os.system('rm ' + save_file)

# Remove some unnecessary files
os.system('rm ' + hits_file)
os.system('rm ' + background_file)
import numpy as np 
import os
import sys
import pdb
import gzip
import scipy.stats



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
    dicti = {}
    for ele in gene_symbol_hits:
    	dicti[ele] = 1
    return dicti

def extract_cardiomyopathy_gene_set(cardiomyopathy_gene_list, gencode_file, version):
	valid_gene_symbols = {}
	f = gzip.open(gencode_file)
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split()
		if line.startswith('#'):
			continue
		line_gene_symbol = data[17].split('"')[1]
		valid_gene_symbols[line_gene_symbol] = 1
	f.close()
	f = open(cardiomyopathy_gene_list)
	head_count = 0
	cardio_genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		geney = data[0]
		if geney not in valid_gene_symbols:
			continue
		if version == 'all':
			cardio_genes[geney] = 1
		elif version == data[1]:
			cardio_genes[geney] = 1
	return cardio_genes


def extract_dynamic_qtl_egenes(egene_file):
	f = open(egene_file)
	egenes = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		egenes[data[1]] = 1
	return egenes

def extract_top_n_dynamic_qtl_egenes(egene_file, num_genes):
	f = open(egene_file)
	egenes = {}
	listy = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		pvalue = float(data[-1])
		ensamble_id = data[1]
		listy.append((ensamble_id, pvalue))
	f.close()
	listy.sort(key=lambda x: x[1])
	for i, tupler in enumerate(listy):
		if i < num_genes:
			egenes[tupler[0]] = 1
	return egenes


def enrichment_analysis(gene_symbol_hits, gene_symbol_background, cardiomyopathy_gene_set, gene_set_name, t):
	aa = 0
	bb = 0
	cc = 0
	dd = 0
	genes = {}
	for gene_id in gene_symbol_hits.keys():
		if gene_id in cardiomyopathy_gene_set:
			genes[gene_id] = 1
			aa = aa + 1
		else:
			bb = bb + 1
	if aa + bb != len(gene_symbol_hits):
		print('ASSUMPTION ERROROR')
	for gene_id in gene_symbol_background.keys():
		if gene_id in gene_symbol_hits:
			continue
		if gene_id in cardiomyopathy_gene_set:
			cc = cc + 1
		else:
			dd = dd + 1
	gene_arr = []
	for gene in genes.keys():
		gene_arr.append(gene)
	gene_arr = np.asarray(gene_arr)
	orat, pvalue = scipy.stats.fisher_exact([[aa,bb],[cc,dd]])
	t.write(gene_set_name + '\t' + str(aa) + '\t' + str(bb) + '\t' + str(cc) + '\t' + str(dd) + '\t' + str(orat) + '\t' + str(pvalue) + '\t' + ','.join(gene_arr) + '\n')
	return t


parameter_string = sys.argv[1]
egene_file = sys.argv[2]
gencode_file = sys.argv[3]
gene_set_enrichment_dir = sys.argv[4]
cardiomyopathy_gene_list = sys.argv[5]
parameter_string = sys.argv[6]
all_association_file = sys.argv[7]

gene_sets = ['all', 'arrhythmogenic_rv_cardiomyopathy', 'dilated_cardiomyopathy_genes', 'hypertrophic_cardiomyopathy_genes', 'lv_non_compaction_cardiomyopathy', 'metabolic_cardiomyopathy_genes']

gene_set_dicti = {}
for gene_set in gene_sets:
	print(gene_set)
	cardiomyopathy_gene_set = extract_cardiomyopathy_gene_set(cardiomyopathy_gene_list, gencode_file, gene_set)
	gene_set_dicti[gene_set] = cardiomyopathy_gene_set

ensamble_hits = extract_dynamic_qtl_egenes(egene_file)
ensamble_background = extract_dynamic_qtl_egenes(all_association_file)
gene_symbol_hits = convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file)
gene_symbol_background = convert_from_ensamble_to_gene_symbol(ensamble_background, gencode_file)
output_file = gene_set_enrichment_dir + parameter_string + '_sig_genes_cardio_gene_set_enrichment.txt'
t = open(output_file, 'w')
t.write('gene_set_version\taa\tbb\tcc\tdd\torat\tpvalue\tgenes\n')
for gene_set in gene_sets:
	t = enrichment_analysis(gene_symbol_hits, gene_symbol_background, gene_set_dicti[gene_set], gene_set, t)
t.close()




num_genes = 200
ensamble_hits = extract_top_n_dynamic_qtl_egenes(egene_file, num_genes)
ensamble_background = extract_dynamic_qtl_egenes(all_association_file)
gene_symbol_hits = convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file)
gene_symbol_background = convert_from_ensamble_to_gene_symbol(ensamble_background, gencode_file)
output_file = gene_set_enrichment_dir + parameter_string + '_top_' + str(num_genes) + '_cardio_gene_set_enrichment.txt'
t = open(output_file, 'w')
t.write('gene_set_version\taa\tbb\tcc\tdd\torat\tpvalue\tgenes\n')
for gene_set in gene_sets:
	t = enrichment_analysis(gene_symbol_hits, gene_symbol_background, gene_set_dicti[gene_set], gene_set, t)
t.close()
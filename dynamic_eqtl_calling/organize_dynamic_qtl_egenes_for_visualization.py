import numpy as np 
import os
import sys
import pdb



# Extract list of of eqtls
def extract_significant_eqtls(significant_qtl_file):
	f = open(significant_qtl_file)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		rs_id = data[0]
		ensamble_id = data[1]
		dicti[ensamble_id + '_' + rs_id] = 1
	f.close()
	dicti['ENSG00000115641_rs11124033'] = 1
	dicti['ENSG00000183873_rs7633988'] = 1
	dicti['ENSG00000183873_rs6599234'] = 1
	dicti['ENSG00000166704_rs8107849'] = 1
	dicti['ENSG00000167173_rs28818910'] = 1
	return dicti




#####################
# Command Line Args
######################
dynamic_eqtl_input_file = sys.argv[1]
significant_qtl_file = sys.argv[2]
significant_egene_file = sys.argv[3]
output_root = sys.argv[4]



output_file = output_root + '_dynamic_qtl_efdr_05_visualization_input.txt'
############
# Extract list of of eqtls
variant_gene_pairs = extract_significant_eqtls(significant_qtl_file)
###########################################
# Filter $dynamic_eqtl_input_file to only those in variant_gene_pairs
###########################################
f = open(dynamic_eqtl_input_file)
t = open(output_file, 'w')
for line in f:
	line = line.rstrip()
	data = line.split()
	rs_id = data[0]
	ensamble_id = data[1]
	if ensamble_id + '_' + rs_id in variant_gene_pairs:
		t.write(line + '\n')
t.close()
f.close()

output_file = output_root + '_dynamic_egene_efdr_05_visualization_input.txt'
############
# Extract list of of eqtls
variant_gene_pairs = extract_significant_eqtls(significant_egene_file)
###########################################
# Filter $dynamic_eqtl_input_file to only those in variant_gene_pairs
###########################################
f = open(dynamic_eqtl_input_file)
t = open(output_file, 'w')
for line in f:
	line = line.rstrip()
	data = line.split()
	rs_id = data[0]
	ensamble_id = data[1]
	if ensamble_id + '_' + rs_id in variant_gene_pairs:
		t.write(line + '\n')
t.close()
f.close()

import numpy as np 
import os
import sys
import pdb
import gzip



def extract_significant_dynamic_qtls(eqtl_file):
	rs_ids = {}
	f = open(eqtl_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		rs_ids[data[0]] = 1
	return rs_ids



def check_for_overlaps(gwas_dir, rs_id_dicti, gwas_threshold, output_file):
	t = open(output_file, 'w')
	t.write('rs_id\thg38_gtex_snp_id\tpvalue\tphenotype\n')
	for file_name in sorted(os.listdir(gwas_dir)):
		if file_name.endswith('.txt') == False:
			continue
		phenotype = file_name.split('.txt')[0]
		print(phenotype)
		gwas_file = gwas_dir + file_name
		f = open(gwas_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				continue
			rs_id = data[0]
			pvalue = data[7]
			if pvalue == 'NA':
				continue
			if float(pvalue) > gwas_threshold:
				continue
			if rs_id not in rs_id_dicti:
				continue
			t.write(rs_id + '\t' + data[1] + '\t' + str(pvalue) + '\t' + phenotype + '\n')
		t.flush()
		f.close()
	t.close()


def extract_significant_per_time_step_qtls(eqtl_file):
	f = open(eqtl_file)
	head_count = 0
	rs_ids = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		rs_ids[data[3]] = 1
	f.close()
	return rs_ids



gtex_gwas_hits_dir = sys.argv[1]
gwas_overlap_output_dir = sys.argv[2]
eqtl_file = sys.argv[3]
parameter_string = sys.argv[4]
gwas_threshold = sys.argv[5]

rs_ids = extract_significant_dynamic_qtls(eqtl_file)
check_for_overlaps(gtex_gwas_hits_dir, rs_ids, float(gwas_threshold), gwas_overlap_output_dir + parameter_string + '_gtex_gwas_' + gwas_threshold + '_overlaps.txt')
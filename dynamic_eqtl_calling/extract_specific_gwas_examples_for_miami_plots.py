import numpy as np 
import os
import sys
import pdb

def mapping_from_hg19_rsid_to_chrom_and_pos(genotype_file):
	f = open(genotype_file)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if line.startswith('#'):
			continue
		rs_id = data[2]
		chrom_num = data[0]
		pos = int(data[1])
		dicti[rs_id] = (chrom_num, pos)
	return dicti

def dynamic_qtl_pvalue_distribution(variant_chrom, variant_pos, output_file, all_hits_file,hg19_rs_id_mapping, ensamble_id):
	t = open(output_file, 'w')
	t.write('rs_id\tvariant_chrom\tvariant_position\tpvalue\n')
	start = 100000000000000000000000
	end = -100000
	f = open(all_hits_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if data[1] != ensamble_id:
			continue
		rs_id = data[0]
		info = hg19_rs_id_mapping[rs_id]
		curr_chrom = info[0]
		curr_pos = int(info[1])
		if curr_pos <= start:
			start = curr_pos
		if curr_pos >= end:
			end = curr_pos
		pvalue = data[-1]
		if curr_chrom != variant_chrom:
			print('ASSUMPTIONE ROROERO')
			pdb.set_trace()
		t.write(rs_id + '\t' + curr_chrom + '\t' + str(curr_pos) + '\t' + pvalue + '\n')
	t.close()
	f.close()
	print(start)
	print(end)
	return start, end

def make_input_bed_file(gwas_file, temporary_hg38_bed_file, variant_chrom):
	t = open(temporary_hg38_bed_file,'w')
	f = open(gwas_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[1] == 'NA':
			continue
		info = data[1].split('_')
		if len(info) == 1:
			pdb.set_trace()
		chromer = info[0]
		pos = info[1]
		if 'chr' + variant_chrom != chromer:
			continue
		t.write(chromer + '\t' + pos + '\t' + str(int(pos) + 1) + '\n')
	t.close()
	f.close()


#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file, output_file, missing_file, liftover_directory):
    stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + 'hg38ToHg19.over.chain.gz ' + output_file + ' ' + missing_file
    os.system(stringer)

#Some jxns were not able to be mapped with liftover (lacked confidence). So first extract those unmapped jxns
#Create dictionary of jxn ids that were not able to mapped via liftover with confidence
def get_unmapped_jxns(temporary_missing_file):
    f = open(temporary_missing_file)
    unmapped_jxns = {}
    for line in f:
        if line.startswith('#'):  # Skip header files
            continue
        line = line.rstrip()
        data = line.split()
        jxn_id = data[0] + ':' + data[1] + ':' + data[2]
        unmapped_jxns[jxn_id] = 1
    return unmapped_jxns

#If raw_leafcutter_file has n+1 lines (therefor n jxns), compute binary vector of length n called pass_filter
#If pass_filter == 1 then that jxn has a mapping in hg19:
def determine_jxns_that_pass_unmappable_filter(cluster_file_hg38, unmapped_jxns):
    pass_filter = []
    count = 0  # for header
    total_jxns = 0  # keep track of number of jxns in the file
    unmapped_jxn_count = 0  # keep track of unmapped jxns
    cluster_counts = {}  # keep track of number of times a cluster is called
    valid_clusters = {}  # only keep clusters that have more than one jxn (after unmapping filter)
    f = open(cluster_file_hg38)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if count == 0:  # skip header
            count = count + 1
            continue
        total_jxns = total_jxns + 1
        jxn_id_info = data[0].split(':')
        cluster_id = jxn_id_info[-1]
        jxn_identifier = jxn_id_info[0] + ':' + jxn_id_info[1] + ':' + jxn_id_info[2]
        if jxn_identifier in unmapped_jxns:  # jxn was unmapped
            pass_filter.append(0)
            unmapped_jxn_count = unmapped_jxn_count + 1
        else:  # jxn mapped
            pass_filter.append(1)
            if cluster_id not in cluster_counts:
                cluster_counts[cluster_id] = 1
            else:
                cluster_counts[cluster_id] = cluster_counts[cluster_id] + 1
    if len(pass_filter) != total_jxns:  # simple check
        print('Length error')
        pdb.set_trace()
    if unmapped_jxn_count != len(unmapped_jxns):  # simple check
        print('lerngther errror')
        pdb.set_trace()
    for cluster_id in cluster_counts.keys():
        if cluster_counts[cluster_id] > 1:  # more than 1 jxn are in this cluster after filter
            valid_clusters[cluster_id] = 1
    return pass_filter, valid_clusters

#Extract converted hg19 coordinates
def get_hg19_coordinates(temporary_hg19_bed_file):
    hg19_coordinates = []
    f = open(temporary_hg19_bed_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        cluster_id = data[0] + ':' + data[1] + ':' + data[2]
        hg19_coordinates.append(cluster_id)
    return hg19_coordinates

def make_hg38_to_hg19_mapping(temporary_hg38_bed_file, temporary_hg19_bed_file):
	f = open(temporary_hg38_bed_file)
	g = open(temporary_hg19_bed_file)
	dicti = {}
	for hg38_line in f:
		hg38_data = hg38_line.rstrip().split()
		hg19_data = g.next().rstrip().split()
		if hg38_data[0] != hg19_data[0]:
			print('assumptione rororo')
		dicti[hg38_data[1]] = hg19_data[1]
	f.close()
	g.close()
	return dicti

def run_liftover_shell(gwas_file, gwas_liftover_root,liftover_directory, variant_chrom, start, end, output_file):
	#First make tempory bed file of our existing hg38 cluster positions. Will be used as input to liftover
	temporary_hg38_bed_file = gwas_liftover_root + 'temp_hg38.bed'
	make_input_bed_file(gwas_file, temporary_hg38_bed_file, variant_chrom)
	temporary_hg19_bed_file = gwas_liftover_root + 'temp_hg19.bed'  # temporary liftover output file
	temporary_missing_file = gwas_liftover_root + 'temp_liftover_missing.bed'  # temporary liftover missing values file
	run_liftover(temporary_hg38_bed_file, temporary_hg19_bed_file, temporary_missing_file, liftover_directory)
    #Some jxns were not able to be mapped with liftover (lacked confidence). So first extract those unmapped jxns
    #Create dictionary where keys are jxn ids that were not able to be mapped via liftover
	unmapped_jxns = get_unmapped_jxns(temporary_missing_file)
	if len(unmapped_jxns) != 0:
		print('ASSSUMPTERION ERROROR')
		pdb.set_trace()
	hg38_to_hg19_mapping = make_hg38_to_hg19_mapping(temporary_hg38_bed_file, temporary_hg19_bed_file)
	os.system('rm ' + temporary_hg38_bed_file)
	os.system('rm ' + temporary_hg19_bed_file)
	os.system('rm ' + temporary_missing_file)

	t = open(output_file, 'w')
	t.write('rs_id\tvariant_chrom\tvariant_position\tpvalue\n')
	f = open(gwas_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[1] == 'NA':
			continue
		info = data[1].split('_')
		chromer = info[0]
		if 'chr' + variant_chrom != chromer:
			continue
		pos = info[1]
		hg19_pos = int(hg38_to_hg19_mapping[pos])
		if hg19_pos >= start and hg19_pos <= end:
			t.write(data[0] + '\t' + variant_chrom + '\t' + str(hg19_pos) + '\t' + data[7] + '\n')
	print('done')
	t.close()
	f.close()

def gwas_pvalue_distribution_neale(gwas_file, variant_chrom, start, end, output_file):
	t = open(output_file, 'w')
	t.write('rs_id\tvariant_chrom\tvariant_position\tpvalue\n')
	f = open(gwas_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		temp_line_chrom = data[0].split(':')[0]
		ref = data[0].split(':')[2]
		alt = data[0].split(':')[3]
		if len(ref) != 1 or len(alt) != 1:
			continue

		if variant_chrom != temp_line_chrom:
			continue
		temp_variant_pos = int(data[0].split(':')[1])
		if temp_variant_pos >= start and temp_variant_pos <= end:
			pvalue = data[-1]
			snp_id = data[0]
			t.write(snp_id + '\t' + variant_chrom + '\t' + str(temp_variant_pos) + '\t' + str(pvalue) + '\n')
	f.close()
	t.close()




def get_pvalue_distributions(rs_id, ensamble_id, variant_chrom, variant_pos, phenotype_name, hg19_rs_id_mapping, all_hits_file, gtex_gwas_hits_dir, gwas_overlap_dir, liftover_directory):
	# First get dynamic qtl pvalue distribution
	output_file = gwas_overlap_dir + rs_id + '_' + ensamble_id + '_' + phenotype_name + '_nearby_dynamic_qtl_pvalues.txt'
	start, end = dynamic_qtl_pvalue_distribution(variant_chrom, variant_pos, output_file, all_hits_file, hg19_rs_id_mapping, ensamble_id)

	gwas_file = gtex_gwas_hits_dir + phenotype_name + '.txt'
	gwas_liftover_root = gwas_overlap_dir + phenotype_name + '_' + variant_chrom + '_hg19_'
	output_file = gwas_overlap_dir + rs_id + '_' + ensamble_id + '_' + phenotype_name + '_nearby_gwas_pvalues.txt'
	run_liftover_shell(gwas_file, gwas_liftover_root,liftover_directory, variant_chrom, start, end, output_file)

def get_pvalue_distributions_from_neil_data(rs_id, ensamble_id, variant_chrom, variant_pos, phenotype_name, hg19_rs_id_mapping, all_hits_file, gwas_file, gwas_overlap_dir):
	# First get dynamic qtl pvalue distribution
	output_file = gwas_overlap_dir + rs_id + '_' + ensamble_id + '_' + phenotype_name + '_nearby_dynamic_qtl_pvalues.txt'
	start, end = dynamic_qtl_pvalue_distribution(variant_chrom, variant_pos, output_file, all_hits_file, hg19_rs_id_mapping, ensamble_id)

	output_file = gwas_overlap_dir + rs_id + '_' + ensamble_id + '_' + phenotype_name + '_nearby_gwas_pvalues.txt'
	gwas_pvalue_distribution_neale(gwas_file, variant_chrom, start, end, output_file)

gtex_gwas_hits_dir = sys.argv[1]
gwas_overlap_dir = sys.argv[2]
all_hits_file = sys.argv[3]
genotype_file = sys.argv[4]
liftover_directory = sys.argv[5]

hg19_rs_id_mapping = mapping_from_hg19_rsid_to_chrom_and_pos(genotype_file)

rs_id = 'rs28818910'
ensamble_id = 'ENSG00000167173'
variant_chrom = '15'
variant_pos = 75440669
phenotype_name = 'UKB_21001_Body_mass_index_BMI'


get_pvalue_distributions(rs_id, ensamble_id, variant_chrom, variant_pos, phenotype_name, hg19_rs_id_mapping, all_hits_file, gtex_gwas_hits_dir, gwas_overlap_dir, liftover_directory)


rs_id = 'rs28818910'
ensamble_id = 'ENSG00000167173'
variant_chrom = '15'
variant_pos = 75440669
phenotype_name = 'Astle_et_al_2016_Red_blood_cell_count'

get_pvalue_distributions(rs_id, ensamble_id, variant_chrom, variant_pos, phenotype_name, hg19_rs_id_mapping, all_hits_file, gtex_gwas_hits_dir, gwas_overlap_dir, liftover_directory)


rs_id = 'rs28818910'
ensamble_id = 'ENSG00000167173'
variant_chrom = '15'
variant_pos = 75440669
phenotype_name = 'UKB_23099_Body_fat_percentage'

get_pvalue_distributions(rs_id, ensamble_id, variant_chrom, variant_pos, phenotype_name, hg19_rs_id_mapping, all_hits_file, gtex_gwas_hits_dir, gwas_overlap_dir, liftover_directory)






#############################
# OLD (RETIRED) SCRIPTS
#############################
'''

rs_id = 'rs28818910'
ensamble_id = 'ENSG00000167173'
variant_chrom = '15'
variant_pos = 75440669
phenotype_name = 'UKB_21001_Body_mass_index_BMI'
gwas_file = gtex_gwas_hits_dir + '21001.assoc.tsv'

get_pvalue_distributions_from_neil_data(rs_id, ensamble_id, variant_chrom, variant_pos, phenotype_name, hg19_rs_id_mapping, all_hits_file, gwas_file, gwas_overlap_dir)


rs_id = 'rs28818910'
ensamble_id = 'ENSG00000167173'
variant_chrom = '15'
variant_pos = 75440669
phenotype_name = 'UKB_23099_Body_fat_percentage'
gwas_file = gtex_gwas_hits_dir + '23099.assoc.tsv'

get_pvalue_distributions_from_neil_data(rs_id, ensamble_id, variant_chrom, variant_pos, phenotype_name, hg19_rs_id_mapping, all_hits_file, gwas_file, gwas_overlap_dir)

'''

import numpy as np 
import os
import sys
import pdb
import gzip




# Get dictionary mapping from cell line to cell line assignment
def get_cell_line_assignments(cell_line_cluster_assignment_file):
	assignment = {}
	f = open(cell_line_cluster_assignment_file)
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		cell_line_name = data[0]
		cluster_assignment = data[1]
		if cluster_assignment == '0':
			cluster_assignment = '2'
		assignment[cell_line_name] = cluster_assignment
	return assignment


# Get mapping from gene symbols to ensamble ids
def get_mapping_from_gene_symbol_to_ensamble_id(gencode_gene_annotation_file,valid_ensamble_ids):
	mapping = {}
	f = gzip.open(gencode_gene_annotation_file)
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split()
		if line.startswith('#'):
			continue
		line_ensamble_id = data[9].split('"')[1].split('.')[0]
		line_gene_symbol = data[17].split('"')[1]
		if line_ensamble_id not in valid_ensamble_ids:
			continue
		mapping[line_gene_symbol] = line_ensamble_id
	return mapping

def get_list_of_used_ensamble_ids(expression_file):
	f = open(expression_file)
	head_count = 0
	used_genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		used_genes[data[0]] = 1
	return used_genes

def get_ensamble_id_cluster_assignments(gene_cluster_assignment_file):
	mapping = {}
	f = open(gene_cluster_assignment_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[0]
		cluster_assignment = data[1]
		mapping[ensamble_id] = str(int(cluster_assignment) + 1)
	return mapping

# Get mapping from each time step to list of column indices that corresponds to cell lines in current cluster at this time step
def get_mapping_from_time_step_to_column_indices(expression_file, cell_line_cluster, cell_line_assignments):
	# First get header line
	f = open(expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			header = data
			head_count = head_count + 1
			continue
		break
	f.close()
	mapping = {}
	# Loop through time steps
	for time_step in range(16):
		valid_indices = []
		for i, ele in enumerate(header):
			if ele == 'Gene_id':
				continue
			cell_line_curr = ele.split('_')[0]
			time_step_curr = ele.split('_')[1]
			if time_step_curr == str(time_step) and cell_line_assignments[cell_line_curr] == cell_line_cluster:
				valid_indices.append(i)
		valid_indices = np.asarray(valid_indices)
		mapping[time_step] = valid_indices
	return mapping

# Extract number of gene cluster
def extract_number_of_gene_clusters(ensamble_assignments):
	assignments = {}
	for ensamble_id in ensamble_assignments.keys():
		assignment = ensamble_assignments[ensamble_id]
		assignments[assignment] = 1
	return len(assignments)

def get_average_expression_for_a_cell_line_cluster(cell_line_cluster, cell_line_assignments, ensamble_assignments, expression_file, output_file):
	# Get mapping from each time step to list of column indices that corresponds to cell lines in current cluster at this time step
	time_steps_to_cell_line_indices = get_mapping_from_time_step_to_column_indices(expression_file, cell_line_cluster, cell_line_assignments)
	# Extract number of gene cluster
	num_gene_clusters = extract_number_of_gene_clusters(ensamble_assignments)
	
	# Open output file and write header
	t = open(output_file, 'w')
	t.write('ensamble_id\tgene_cluster')
	for time_step in range(16):
		t.write('\t' + 'average_expression_' + str(time_step))
	t.write('\n')

	# Loop through gene clusters
	for gene_cluster_num in range(1, num_gene_clusters+1):
		gene_cluster_name = str(gene_cluster_num)
		print(gene_cluster_num)

		# Loop through expression file for this gene cluster
		f = open(expression_file)
		head_count = 0
		for line in f:
			line = line.strip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				continue
			ensamble_id = data[0]
			if ensamble_id not in ensamble_assignments:
				continue
			if ensamble_assignments[ensamble_id] != gene_cluster_name:
				continue
			t.write(ensamble_id + '\t' + ensamble_assignments[ensamble_id])
			for time_step in range(16):
				average_expression = np.mean(np.asarray(data)[time_steps_to_cell_line_indices[time_step]].astype(float))
				t.write('\t' + str(average_expression))
			t.write('\n')
		f.close()
	t.close()

preprocess_total_expression_dir = sys.argv[1]
mixutre_hmm_cell_line_grouping_dir = sys.argv[2]

# File containing all gene expression measurements
expression_file = preprocess_total_expression_dir + 'quantile_normalized.txt'
# File containing cell line cluster assignments
cell_line_cluster_assignment_file = mixutre_hmm_cell_line_grouping_dir + 'mixsvgp_K2_L100_1_28542829_assignments'
# File containing gene cluster assignments
gene_cluster_assignment_file = mixutre_hmm_cell_line_grouping_dir + 'mixsvgp_K2_L20_0_29526888_gene_assignments'


# Get dictionary mapping from cell line to cell line assignment
cell_line_assignments = get_cell_line_assignments(cell_line_cluster_assignment_file)


# Get dictionary mapping from ensamble id to gene cluster assignment
ensamble_assignments = get_ensamble_id_cluster_assignments(gene_cluster_assignment_file)

######################
# Should add check here once karl gets revised gene assignments
######################

cell_line_cluster = '1'
# Output file for mean expression matrix belonging to cell line cluster 1
output_file = preprocess_total_expression_dir + 'cell_line_cluster' + cell_line_cluster + '_average_expression.txt'
get_average_expression_for_a_cell_line_cluster(cell_line_cluster, cell_line_assignments, ensamble_assignments, expression_file, output_file)


cell_line_cluster = '2'
# Output file for mean expression matrix belonging to cell line cluster 1
output_file = preprocess_total_expression_dir + 'cell_line_cluster' + cell_line_cluster + '_average_expression.txt'
get_average_expression_for_a_cell_line_cluster(cell_line_cluster, cell_line_assignments, ensamble_assignments, expression_file, output_file)

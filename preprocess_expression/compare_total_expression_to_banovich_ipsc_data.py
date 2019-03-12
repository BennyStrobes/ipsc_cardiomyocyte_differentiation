import numpy as np 
import os
import sys
import pdb
import scipy.stats
from sklearn import linear_model


# Extract dictionary list of overlapping samples between two data sets
def extract_list_of_overlapping_samples_between_two_data_sets(banovich_read_counts_file, time_series_read_counts_file, day, data_type):
	time_series_samples = {}
	banovich_samples = {}
	union_samples = {}
	# Extract samples from time series file
	f = open(time_series_read_counts_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for sample_id in data[1:]:
				cell_line_id = sample_id.split('_')[0]
				time_step = sample_id.split('_')[1]
				if time_step == str(day):
					time_series_samples[cell_line_id] = 1
			continue
	f.close()
	# Extract samples from time series file
	f = open(banovich_read_counts_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for sample_id in data[1:]:
				if data_type == 'ipsc_cm':
					sample_id = sample_id.split('A')[1]
				banovich_samples[sample_id] = 1
			continue
	f.close()
	# Take union of two dictionaries
	for sampler in time_series_samples.keys():
		if sampler in banovich_samples:
			union_samples[sampler] = 1
	return union_samples

def get_good_genes_standard_pca(gene_weighs_file):
	good_genes = {}
	f = open(gene_weighs_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[0]
		pc1 = float(data[1])
		pc2 = float(data[2])
		pc3 = float(data[3])
		pc4 = float(data[4])
		pc5 = float(data[5])
		pc6 = float(data[6])
		pc7 = float(data[7])
		threshold = .005
		if abs(pc1) < threshold and abs(pc2) < threshold and abs(pc3) < threshold and abs(pc4) < threshold and abs(pc5) < threshold:
			good_genes[ensamble_id] = 1
	f.close()
	print(len(good_genes))
	return good_genes

def get_good_genes_based_on_cell_line_pca(day, file_name):
	good_genes = {}
	f = open(file_name)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_namer = data[0]
		ensamble_id = gene_namer.split('_')[0]
		line_day = gene_namer.split('_')[1]
		if line_day != day:
			continue
		pc1 = float(data[1])
		pc2 = float(data[2])
		pc3 = float(data[3])
		pc4 = float(data[4])
		pc5 = float(data[5])
		pc6 = float(data[6])
		pc7 = float(data[7])
		threshold = .002
		if abs(pc1) < threshold and abs(pc2) < threshold and abs(pc3) < threshold and abs(pc4) < threshold and abs(pc5) < threshold:
			good_genes[ensamble_id] = 1
	f.close()
	print(len(good_genes))



	return good_genes

# Compute the correlation
def take_correlation_between_two_data_types(banovich_sample, time_series_sample, banovich_read_counts_file, time_series_read_counts_file, day, data_type, good_genes):
	# good_genes = get_good_genes_standard_pca('/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/covariates/principal_components_10_gene_weights.txt')
	# good_genes = get_good_genes_based_on_cell_line_pca(str(day), '/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/covariates/cell_line_ignore_missing_principal_components_10_gene_weights.txt')
	# Create mapping from gene id to raw read counts
	banovich_genes = {}
	time_series_genes = {}
	# Fill in this mapping for time series file
	f = open(time_series_read_counts_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for i, sample_id in enumerate(data[1:]):
				cell_line_id = sample_id.split('_')[0]
				time_step = sample_id.split('_')[1]
				if time_step == str(day) and cell_line_id == time_series_sample:
					column_num = i + 1
			continue
		ensamble_id = data[0]
		read_count = float(data[column_num])
		if ensamble_id in time_series_genes:
			print('assumption error!')
			pdb.set_trace()
		time_series_genes[ensamble_id] = read_count
	f.close()
	# Fill in this mapping for banovich data
	f = open(banovich_read_counts_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for i, sample_id in enumerate(data[1:]):
				if data_type == 'ipsc_cm':
					sample_id = sample_id.split('A')[1]
				if sample_id == banovich_sample:
					column_num = i + 1
			continue
		ensamble_id = data[0].split('.')[0]
		read_count = float(data[column_num])
		if ensamble_id in banovich_genes:
			print('assumption error!!')
			pdb.set_trace()
		banovich_genes[ensamble_id] = read_count
	f.close()
	# Now put read counts into array where each element corresponds to the same gene (only if that gene is presenent in both data sets)
	time_series_arr = []
	banovich_arr = []
	for key in banovich_genes.keys():
		if key in time_series_genes:
			if key not in good_genes:
				continue
			banovich_arr.append(banovich_genes[key])
			time_series_arr.append(time_series_genes[key])
	rho, pvalue = scipy.stats.pearsonr(np.log2(np.asarray(banovich_arr)+1), np.log2(np.asarray(time_series_arr)+1))
	return rho

# Compute the correlation
def take_correlation_between_two_data_types_regress_out_version(banovich_sample, time_series_sample, banovich_read_counts_file, time_series_read_counts_file):
	# good_genes = get_good_genes_standard_pca('/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/covariates/principal_components_10_gene_weights.txt')
	# good_genes = get_good_genes_based_on_cell_line_pca(str(day), '/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/covariates/cell_line_ignore_missing_principal_components_10_gene_weights.txt')
	# Create mapping from gene id to raw read counts
	banovich_genes = {}
	time_series_genes = {}
	# Fill in this mapping for time series file
	f = open(time_series_read_counts_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for i, cell_line_id in enumerate(data[1:]):
				if cell_line_id == time_series_sample:
					column_num = i + 1
			continue
		ensamble_id = data[0]
		read_count = float(data[column_num])
		if ensamble_id in time_series_genes:
			print('assumption error!')
			pdb.set_trace()
		time_series_genes[ensamble_id] = read_count
	f.close()
	# Fill in this mapping for banovich data
	f = open(banovich_read_counts_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for i, sample_id in enumerate(data[1:]):
				if sample_id == banovich_sample:
					column_num = i + 1
			continue
		ensamble_id = data[0]
		read_count = float(data[column_num])
		if ensamble_id in banovich_genes:
			print('assumption error!!')
			pdb.set_trace()
		banovich_genes[ensamble_id] = read_count
	f.close()
	# Now put read counts into array where each element corresponds to the same gene (only if that gene is presenent in both data sets)
	time_series_arr = []
	banovich_arr = []
	for key in banovich_genes.keys():
		if key in time_series_genes:
			banovich_arr.append(banovich_genes[key])
			time_series_arr.append(time_series_genes[key])
	rho, pvalue = scipy.stats.pearsonr(np.asarray(banovich_arr), np.asarray(time_series_arr))
	return rho

# Compute the correlation
def take_correlation_between_two_data_types_all_genes(banovich_sample, time_series_sample, banovich_read_counts_file, time_series_read_counts_file, day, data_type):
	# good_genes = get_good_genes_standard_pca('/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/covariates/principal_components_10_gene_weights.txt')
	# good_genes = get_good_genes_based_on_cell_line_pca(str(day), '/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/covariates/cell_line_ignore_missing_principal_components_10_gene_weights.txt')
	# Create mapping from gene id to raw read counts
	banovich_genes = {}
	time_series_genes = {}
	# Fill in this mapping for time series file
	f = open(time_series_read_counts_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for i, sample_id in enumerate(data[1:]):
				cell_line_id = sample_id.split('_')[0]
				time_step = sample_id.split('_')[1]
				if time_step == str(day) and cell_line_id == time_series_sample:
					column_num = i + 1
			continue
		ensamble_id = data[0]
		read_count = float(data[column_num])
		if ensamble_id in time_series_genes:
			print('assumption error!')
			pdb.set_trace()
		time_series_genes[ensamble_id] = read_count
	f.close()
	# Fill in this mapping for banovich data
	f = open(banovich_read_counts_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for i, sample_id in enumerate(data[1:]):
				if data_type == 'ipsc_cm':
					sample_id = sample_id.split('A')[1]
				if sample_id == banovich_sample:
					column_num = i + 1
			continue
		ensamble_id = data[0].split('.')[0]
		read_count = float(data[column_num])
		if ensamble_id in banovich_genes:
			print('assumption error!!')
			pdb.set_trace()
		banovich_genes[ensamble_id] = read_count
	f.close()
	# Now put read counts into array where each element corresponds to the same gene (only if that gene is presenent in both data sets)
	time_series_arr = []
	banovich_arr = []
	for key in banovich_genes.keys():
		if key in time_series_genes:
			banovich_arr.append(banovich_genes[key])
			time_series_arr.append(time_series_genes[key])
	rho, pvalue = scipy.stats.pearsonr(np.log2(np.asarray(banovich_arr)+1), np.log2(np.asarray(time_series_arr)+1))
	return rho

def standardize_rows(counts, genes):
	row, col = counts.shape
	good_rows = []
	new_counts = np.zeros((row,col))
	for row_num in range(row):
		new_counts[row_num,:] = (counts[row_num,:] - np.mean(counts[row_num,:]))/np.std(counts[row_num,:])
		if np.std(counts[row_num,:]) != 0:
			good_rows.append(row_num)
	return new_counts[good_rows,:], genes[good_rows]

def standardize_rows_v2(counts, genes):
	row, col = counts.shape
	good_rows = []
	new_counts = np.zeros((row,col))
	for row_num in range(row):
		new_counts[row_num,:] = (counts[row_num,:] - np.mean(counts[row_num,:]))/np.std(counts[row_num,:])
		if np.std(counts[row_num,:]) != 0:
			good_rows.append(row_num)
	return new_counts[good_rows,:], counts[good_rows,:], genes[good_rows]

def center_rows(counts, genes):
	row, col = counts.shape
	good_rows = []
	new_counts = np.zeros((row,col))
	for row_num in range(row):
		new_counts[row_num,:] = counts[row_num,:] - np.mean(counts[row_num,:])
		if np.std(counts[row_num,:]) != 0:
			good_rows.append(row_num)
	return new_counts[good_rows,:], counts[good_rows,:], genes[good_rows]

def get_good_genes_via_time_independent_pca(u_mat, gene_arr, num_pc, threshold):
	good_genes = {}
	for index in range(len(gene_arr)):
		gene_id = gene_arr[index]
		valid = True 
		for pc_num in range(num_pc):
			if abs(u_mat[index, pc_num]) >= threshold:
				valid = False
		if valid == True:
			good_genes[gene_id] = 1
	return good_genes

def get_good_genes_via_time_independent_pca_v2(u_mat, gene_arr, num_pc, num_genes):
	good_genes = {}
	loadings = np.amax(np.abs(u_mat[:,:num_pc]),axis=1)

	tupler_list = []
	for i, loading in enumerate(loadings):
		tupler_list.append((gene_arr[i], loading))

	
	tupler_list.sort(key=lambda x:x[1])

	for i, ele in enumerate(tupler_list):
		if i < num_genes:
			good_genes[ele[0]] = 1

	return good_genes


def get_good_genes(banovich_read_counts_file, time_series_read_counts_file, day, data_type, num_genes, num_pc_banovich, num_pc_time_series):
	# Get banovich read count matrix
	banovich_count_data = np.loadtxt(banovich_read_counts_file, dtype=str,delimiter='\t')
	banovich_genes = banovich_count_data[1:,0]
	banovich_counts = np.log2(banovich_count_data[1:,1:].astype(float) + 1)
	if data_type == 'ipsc':
		for i, ele in enumerate(banovich_genes):
			banovich_genes[i] = ele.split('.')[0]
	# Get time series read count matrix
	time_series_count_data = np.loadtxt(time_series_read_counts_file, dtype=str)
	time_series_genes = time_series_count_data[1:,0]
	time_series_counts = np.log2(time_series_count_data[1:,1:].astype(float) + 1)
	time_series_labels = time_series_count_data[0,1:]
	good_labels = []
	for i,ele in enumerate(time_series_labels):
		if ele.split('_')[1] == day:
			good_labels.append(i)
	time_series_counts = time_series_counts[:, good_labels]

	standardized_banovich_counts, banovich_genes = standardize_rows(banovich_counts, banovich_genes)
	standardized_time_series_counts,time_series_genes  = standardize_rows(time_series_counts,time_series_genes)

	u_banovich, s_banovich, v_banovich = np.linalg.svd(standardized_banovich_counts)
	u_time, s_time, v_time = np.linalg.svd(standardized_time_series_counts)
	
	#valid_banovich_genes = get_good_genes_via_time_independent_pca(u_banovich, banovich_genes, 3, .008)
	#valid_time_series_genes = get_good_genes_via_time_independent_pca(u_time, time_series_genes, 3, .008)

	valid_banovich_genes = get_good_genes_via_time_independent_pca_v2(u_banovich, banovich_genes, num_pc_banovich, num_genes)
	valid_time_series_genes = get_good_genes_via_time_independent_pca_v2(u_time, time_series_genes, num_pc_time_series, num_genes)


	overlapping_genes = {}
	for key in valid_banovich_genes.keys():
		if key in valid_time_series_genes:
			overlapping_genes[key] = 1
	print(len(overlapping_genes))
	return overlapping_genes





def correlation_between_time_series_and_banovich_samples(banovich_read_counts_file, time_series_read_counts_file, day, data_type, num_genes, num_pc_banovich, num_pc_time_series):
	# Extract dictionary list of overlapping samples between two data sets
	overlapping_samples = extract_list_of_overlapping_samples_between_two_data_sets(banovich_read_counts_file, time_series_read_counts_file, day, data_type)
	ordered_samples = sorted(overlapping_samples.keys())

	# Initialize correlation matrix
	num_samples = len(overlapping_samples)
	correlation_mat = np.zeros((num_samples, num_samples))


	good_genes = get_good_genes(banovich_read_counts_file, time_series_read_counts_file, str(day), data_type, num_genes, num_pc_banovich, num_pc_time_series)
	# Loop through all pairs of samples and take the correlation between the two
	matched = []
	unmatched = []
	for i in range(num_samples):
		for j in range(num_samples):
			time_series_sample = ordered_samples[i]
			banovich_sample = ordered_samples[j]
			# Compute the correlation
			correlation_mat[i,j] = take_correlation_between_two_data_types(banovich_sample, time_series_sample, banovich_read_counts_file, time_series_read_counts_file, day, data_type, good_genes)
			if i == j:
				matched.append(correlation_mat[i,j])
			else:
				unmatched.append(correlation_mat[i,j])
	print(scipy.stats.ranksums(matched,unmatched))
	output_mat = np.zeros((num_samples + 1, num_samples + 1)).astype(str)
	output_mat[1:,1:] = correlation_mat.astype(str)
	output_mat[0,0] = 'Cell_line'
	output_mat[0,1:] = ordered_samples
	output_mat[1:,0] = ordered_samples
	return output_mat

def regress_out(expr, loadings):
	num_genes, num_samples = expr.shape
	residual_expr = np.zeros((num_genes, num_samples))
	for gene_num in range(num_genes):
		y = np.transpose(np.asmatrix(expr[gene_num,:]))
		model = linear_model.LinearRegression(fit_intercept=True)
		modelfit = model.fit(loadings, y)
		beta = modelfit.coef_
		for i,val in enumerate(beta[0]):
			y = y - np.transpose(np.asmatrix(val*loadings[:,i]))
		residual_expr[gene_num,:] = np.squeeze(np.asarray(y))
	return residual_expr


def regress_out_pcs_time_series(time_series_read_counts_file, day, num_pc_time_series, time_series_regress_out_counts_file):
	# Get time series read count matrix
	time_series_count_data = np.loadtxt(time_series_read_counts_file, dtype=str)
	time_series_genes = time_series_count_data[1:,0]
	time_series_counts = np.log2(time_series_count_data[1:,1:].astype(float) + 1)
	time_series_labels = time_series_count_data[0,1:]
	good_labels = []
	cell_lines = []
	for i,ele in enumerate(time_series_labels):
		if ele.split('_')[1] == day:
			good_labels.append(i)
			cell_lines.append(ele.split('_')[0])
	time_series_counts = time_series_counts[:, good_labels]

	centered_time_series_counts, time_series_genes  = standardize_rows(time_series_counts,time_series_genes)

	u_time, s_time, v_time = np.linalg.svd(centered_time_series_counts)

	loadings = np.transpose(v_time)[:,:num_pc_time_series]
	
	residual_counts = regress_out(centered_time_series_counts, loadings)

	output_mat = np.zeros((len(time_series_genes) + 1, len(cell_lines) + 1)).astype(str)
	output_mat[1:,1:] = residual_counts.astype(str)
	output_mat[0,0] = 'Gene_id'
	output_mat[0,1:] = cell_lines
	output_mat[1:,0] = time_series_genes
	np.savetxt(time_series_regress_out_counts_file, output_mat, fmt="%s",delimiter='\t')
	return



def regress_out_pcs_banovich(banovich_read_counts_file, data_type, num_pc_banovich, banovich_regress_out_counts_file):
	# Get banovich read count matrix
	banovich_count_data = np.loadtxt(banovich_read_counts_file, dtype=str,delimiter='\t')
	banovich_genes = banovich_count_data[1:,0]
	banovich_counts = np.log2(banovich_count_data[1:,1:].astype(float) + 1)
	banovich_samples = banovich_count_data[0,1:]
	if data_type == 'ipsc':
		for i, ele in enumerate(banovich_genes):
			banovich_genes[i] = ele.split('.')[0]
	samples = banovich_count_data[0,1:]
	if data_type == 'ipsc_cm':
		for i, ele in enumerate(samples):
			samples[i] = ele.split('A')[1]

	centered_banovich_counts, banovich_genes = standardize_rows(banovich_counts, banovich_genes)

	if num_pc_banovich > 0:
		u_banovich, s_banovich, v_banovich = np.linalg.svd(centered_banovich_counts)

		loadings = np.transpose(v_banovich)[:,:num_pc_banovich]
	
		residual_counts = regress_out(centered_banovich_counts, loadings)
	else:
		residual_counts = centered_banovich_counts
	output_mat = np.zeros((len(banovich_genes)+1, banovich_count_data.shape[1])).astype(str)
	output_mat[1:,1:] = residual_counts.astype(str)
	output_mat[0,0] = 'Gene_id'
	output_mat[0,1:] = samples
	output_mat[1:,0] = banovich_genes
	np.savetxt(banovich_regress_out_counts_file, output_mat, fmt="%s",delimiter='\t')
	return

def correlation_between_time_series_and_banovich_samples_regress_out(banovich_read_counts_file, time_series_read_counts_file, day, data_type, num_pc_banovich, num_pc_time_series, output_dir):
	# Extract dictionary list of overlapping samples between two data sets
	overlapping_samples = extract_list_of_overlapping_samples_between_two_data_sets(banovich_read_counts_file, time_series_read_counts_file, day, data_type)
	ordered_samples = sorted(overlapping_samples.keys())

	# Initialize correlation matrix
	num_samples = len(overlapping_samples)
	correlation_mat = np.zeros((num_samples, num_samples))

	time_series_regress_out_counts_file = output_dir + 'time_series_day_' + str(day) + '_' + str(num_pc_time_series) + '_pc_corrected_expr.txt'
	regress_out_pcs_time_series(time_series_read_counts_file, str(day), num_pc_time_series, time_series_regress_out_counts_file)

	banovich_regress_out_counts_file = output_dir + 'banovich_' + data_type + '_' + str(num_pc_banovich) + '_pc_corrected_expr.txt'
	regress_out_pcs_banovich(banovich_read_counts_file, data_type, num_pc_banovich, banovich_regress_out_counts_file)
	# Loop through all pairs of samples and take the correlation between the two
	matched = []
	unmatched = []
	for i in range(num_samples):
		for j in range(num_samples):
			time_series_sample = ordered_samples[i]
			banovich_sample = ordered_samples[j]
			# Compute the correlation
			correlation_mat[i,j] = take_correlation_between_two_data_types_regress_out_version(banovich_sample, time_series_sample, banovich_regress_out_counts_file, time_series_regress_out_counts_file)
			if i == j:
				matched.append(correlation_mat[i,j])
			else:
				unmatched.append(correlation_mat[i,j])
	print(scipy.stats.ranksums(matched,unmatched))
	output_mat = np.zeros((num_samples + 1, num_samples + 1)).astype(str)
	output_mat[1:,1:] = correlation_mat.astype(str)
	output_mat[0,0] = 'Cell_line'
	output_mat[0,1:] = ordered_samples
	output_mat[1:,0] = ordered_samples
	return output_mat


def correlation_between_time_series_and_banovich_samples_all_genes(banovich_read_counts_file, time_series_read_counts_file, day, data_type):
	# Extract dictionary list of overlapping samples between two data sets
	overlapping_samples = extract_list_of_overlapping_samples_between_two_data_sets(banovich_read_counts_file, time_series_read_counts_file, day, data_type)
	ordered_samples = sorted(overlapping_samples.keys())

	# Initialize correlation matrix
	num_samples = len(overlapping_samples)
	correlation_mat = np.zeros((num_samples, num_samples))
	# Loop through all pairs of samples and take the correlation between the two
	matched = []
	unmatched = []
	for i in range(num_samples):
		for j in range(num_samples):
			time_series_sample = ordered_samples[i]
			banovich_sample = ordered_samples[j]
			# Compute the correlation
			correlation_mat[i,j] = take_correlation_between_two_data_types_all_genes(banovich_sample, time_series_sample, banovich_read_counts_file, time_series_read_counts_file, day, data_type)
			if i == j:
				matched.append(correlation_mat[i,j])
			else:
				unmatched.append(correlation_mat[i,j])
	print(scipy.stats.ranksums(matched,unmatched))
	output_mat = np.zeros((num_samples + 1, num_samples + 1)).astype(str)
	output_mat[1:,1:] = correlation_mat.astype(str)
	output_mat[0,0] = 'Cell_line'
	output_mat[0,1:] = ordered_samples
	output_mat[1:,0] = ordered_samples
	return output_mat



######################
# Command line args
######################
ipsc_banovich_read_counts_file = sys.argv[1]
ipsc_cm_banovich_read_counts_file = sys.argv[2]
preprocess_total_expression_dir = sys.argv[3]
banovich_ipsc_comparison_dir = sys.argv[4]


time_series_read_counts_file = preprocess_total_expression_dir + 'raw_counts.txt'




#######################
# Part 1: Take correlation between day 0 time step samples and  $ipsc_banovich_read_counts_file samples after regressing out PCs
num_pc_banovich = 10
num_pc_time_series = 3
output_file = banovich_ipsc_comparison_dir + 'day_0_banovich_ipsc_comparison_regress_out_' + str(num_pc_banovich) + '_' + str(num_pc_time_series) + '_pcs.txt'
correlation_mat = correlation_between_time_series_and_banovich_samples_regress_out(ipsc_banovich_read_counts_file, time_series_read_counts_file, 0, 'ipsc', num_pc_banovich, num_pc_time_series, banovich_ipsc_comparison_dir)
np.savetxt(output_file, correlation_mat, fmt="%s",delimiter='\t')


#######################
# Part 2: Take correlation between day 0 time step samples and  $ipsc_cm_banovich_read_counts_file samples
num_pc_banovich = 3
num_pc_time_series = 3
output_file = banovich_ipsc_comparison_dir + 'day_15_banovich_ipsc_comparison_regress_out_' + str(num_pc_banovich) + '_' + str(num_pc_time_series) + '_pcs.txt'
correlation_mat = correlation_between_time_series_and_banovich_samples_regress_out(ipsc_cm_banovich_read_counts_file, time_series_read_counts_file, 15, 'ipsc_cm', num_pc_banovich, num_pc_time_series, banovich_ipsc_comparison_dir)
np.savetxt(output_file, correlation_mat, fmt="%s",delimiter='\t')






#######################
# Part 1: Take correlation between day 0 time step samples and  $ipsc_banovich_read_counts_file samples
num_genes = 1000
num_pc_banovich = 3
num_pc_time_series=3
output_file = banovich_ipsc_comparison_dir + 'day_0_banovich_ipsc_comparison_' + str(num_pc_banovich) + '_' + str(num_pc_time_series) + '_' + str(num_genes) + '.txt'
correlation_mat = correlation_between_time_series_and_banovich_samples(ipsc_banovich_read_counts_file, time_series_read_counts_file, 0, 'ipsc', num_genes, num_pc_banovich, num_pc_time_series)
np.savetxt(output_file, correlation_mat, fmt="%s",delimiter='\t')
#######################
# Part 2: Take correlation between day 0 time step samples and  $ipsc_cm_banovich_read_counts_file samples
num_genes = 1000
num_pc_banovich = 3
num_pc_time_series=3
output_file = banovich_ipsc_comparison_dir + 'day_15_banovich_ipsc_comparison_' + str(num_pc_banovich) + '_' + str(num_pc_time_series) + '_' + str(num_genes) + '.txt'
correlation_mat = correlation_between_time_series_and_banovich_samples(ipsc_cm_banovich_read_counts_file, time_series_read_counts_file, 15, 'ipsc_cm', num_genes, num_pc_banovich, num_pc_time_series)
np.savetxt(output_file, correlation_mat, fmt="%s",delimiter='\t')




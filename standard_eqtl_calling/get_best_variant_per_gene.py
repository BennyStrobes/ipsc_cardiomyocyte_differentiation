import numpy as np 
import os
import sys
import scipy.stats
import pdb


def extract_all_significant_variant_gene_pairs(qtls, file_name):
    f = open(file_name)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        gene_id = data[1]
        rs_id = data[3]
        test_id = gene_id + '_' + rs_id
        qtls[test_id] = np.zeros(16) + 100
    return qtls


def extract_variant_gene_pairs(cht_output_dir, parameter_string, fdr):
    qtls = {}
    for time_step in range(16):
        file_name = cht_output_dir + 'cht_results_' + parameter_string + '_time_' + str(time_step) + '_efdr_thresh_' + str(fdr) + '_significant.txt'
        qtls = extract_all_significant_variant_gene_pairs(qtls, file_name)

    qtls = fill_in_pvalues(cht_output_dir, parameter_string, qtls)
    return qtls



def fill_in_pvalues(cht_output_dir, parameter_string, qtls):
    for time_step in range(16):
        f = open(cht_output_dir + 'cht_results_' + parameter_string + '_time_' + str(time_step) + '_eqtl_results.txt')
        head_count = 0
        for line in f:
            line = line.rstrip()
            data = line.split()
            if head_count == 0:
                head_count = head_count + 1
                continue
            gene_id = data[1]
            rs_id = data[3]
            test_id = gene_id + '_' + rs_id
            if test_id in qtls:
                pvalue = float(data[8])
                qtls[test_id][time_step] = pvalue
        f.close()
    return qtls





def print_qtls(qtls,output_file):
    t = open(output_file,'w')
    for test_id in qtls.keys():
        pvalues = qtls[test_id]
        pvalue_string = np.asarray(pvalues).astype(str)
        t.write(test_id + '\t' + '\t'.join(pvalue_string) + '\n')
    t.close()

def get_egenes_through_geometric_mean(qtls):
    gene_to_info = {}
    for test_name in qtls.keys():
        gene_id = test_name.split('_')[0]
        variant_id = test_name.split('_')[1]
        pvalue_vector = qtls[test_name]
        # Simple check to make sure all time step's pvalues are filled in
        if len(np.where(pvalue_vector == 100)[0]) > 0:
            print('FATAL ERROR')
            pdb.set_trace()
        # Compute geometric mean
        # Add '.0000000000000000000000001' gmean takes log. And otherwise we get warning for taking log of 0. This does not affect geometric mean
        geometric_mean_pvalue = scipy.stats.mstats.gmean(pvalue_vector + .0000000000000000000000001)

        if gene_id not in gene_to_info: # If gene has never been seen before
            gene_to_info[gene_id] = (variant_id,geometric_mean_pvalue)
        elif gene_id in gene_to_info: # Gene has been seen before
            old_tuple = gene_to_info[gene_id]
            old_variant_id = old_tuple[0]
            old_geometric_pvalue = old_tuple[1]
            if geometric_mean_pvalue < old_geometric_pvalue: # Current geometric mean pvalue is more significant than previous best
                gene_to_info[gene_id] = (variant_id,geometric_mean_pvalue)
            else:  # Prevoius best geometric mean pvalue is more significant than current geometric mean pvalue
                gene_to_info[gene_id] = (old_variant_id, old_geometric_pvalue)

    # Now convert gene_to_info data structure to one where keys are #gene_id + '_' + $variant_id
    # This data structure will be easier to work with in future
    egenes = {}
    for egene in gene_to_info.keys():
        test_name = egene + '_' + gene_to_info[egene][0]
        egenes[test_name] = np.zeros(16) + 100
    return egenes

def print_egene_name_output_file(egenes, egene_name_output_file):
    t = open(egene_name_output_file, 'w')
    t.write('gene_id\tvariant_id\n')
    for test_id in egenes.keys():
        gene_id = test_id.split('_')[0]
        variant_id = test_id.split('_')[1]
        t.write(gene_id + '\t' + variant_id + '\n')
    t.close()

def print_summary_statistic(cht_output_dir, parameter_string, egenes_cp, output_file,column_number):
    for time_step in range(16):
        f = open(cht_output_dir + 'cht_results_' + parameter_string + '_time_' + str(time_step) + '_eqtl_results.txt')
        head_count = 0
        for line in f:
            line = line.rstrip()
            data = line.split()
            if head_count == 0:
                head_count = head_count + 1
                continue
            gene_id = data[1]
            rs_id = data[3]
            test_id = gene_id + '_' + rs_id
            if test_id in egenes:
                statistic = float(data[column_number])
                egenes_cp[test_id][time_step] = statistic
        f.close()
    t = open(output_file, 'w')
    for test_id in sorted(egenes_cp.keys()):
        statistic_vector = egenes_cp[test_id]
        t.write(test_id + '\t' + '\t'.join(statistic_vector.astype(str)) + '\n')
    t.close()


# Create list summarizing number of time steps each variant gene pair is significant in
def eqtl_sharing_result(eqtl_sharing_output, cht_output_dir, parameter_string, fdr):
    t = open(eqtl_sharing_output, 'w')
    hits = {}
    for time_step in range(16):
        file_name = cht_output_dir + 'cht_results_' + parameter_string + '_time_' + str(time_step) + '_efdr_thresh_' + str(fdr) + '_significant.txt'
        f = open(file_name)
        head_count = 0
        for line in f:
            line = line.rstrip()
            data = line.split()
            if head_count == 0:
                head_count = head_count + 1
                continue
            gene = data[1]
            rs_id = data[3]
            test_name = gene + '_' + rs_id
            if test_name not in hits:
                hits[test_name] = []
            hits[test_name].append(time_step)
    f.close()

    for time_step in range(17):
        for hit in hits.keys():
            if len(hits[hit]) == time_step:
                t.write(hit + '\t' + ','.join(np.asarray(hits[hit]).astype(str)) + '\t' + str(time_step) + '\n')
    t.close()


parameter_string = sys.argv[1]
cht_output_dir = sys.argv[2]
num_pc = sys.argv[3]
fdr = sys.argv[4]


parameter_string = parameter_string + '_num_pc_' + num_pc



# Create list summarizing number of time steps each variant gene pair is significant in
eqtl_sharing_output = cht_output_dir + parameter_string +'_fdr_' + str(fdr) + '_eqtl_sharing.txt'
eqtl_sharing_result(eqtl_sharing_output, cht_output_dir, parameter_string, fdr)


# Extract dictionary called qtls with:
## 1. keys that are all $gene_id + '_' + $variant_id (test name) pairs that are genome wide significant in 1 time step
## 2. values are an array of length number of time steps. Where each element in the arrray is the pvalue of that test in that time step
qtls = extract_variant_gene_pairs(cht_output_dir, parameter_string, fdr)

# Create dictionary called egenes with keys:
## 1. that are all $gene_id + '_' + $variant_id (test name) pairs such that the variant has the smallest geometric mean pvalue of all other variants
egenes = get_egenes_through_geometric_mean(qtls)


# Now print results to files
egene_name_output_file = cht_output_dir + 'best_variant_per_egene_' + parameter_string +'_fdr_' + str(fdr) + '_test_names.txt'
print_egene_name_output_file(egenes, egene_name_output_file)

# Now print summary statistics for best variant per gene tests to output file 
# 1 file per summary statistic
# 1 row per test
# Each row has 16 columns (1 per time step)
egene_alpha_output_file = cht_output_dir + 'best_variant_per_egene_' + parameter_string +'_fdr_' + str(fdr) + '_alpha.txt'
print_summary_statistic(cht_output_dir, parameter_string, egenes, egene_alpha_output_file, 6)

egene_beta_output_file = cht_output_dir + 'best_variant_per_egene_' + parameter_string +'_fdr_' + str(fdr) + '_beta.txt'
print_summary_statistic(cht_output_dir, parameter_string, egenes, egene_beta_output_file, 7)

egene_pvalue_output_file = cht_output_dir + 'best_variant_per_egene_' + parameter_string +'_fdr_' + str(fdr) + '_pvalues.txt'
print_summary_statistic(cht_output_dir, parameter_string, egenes, egene_pvalue_output_file, 8)


















######################################################
# OLD STUFF. NO LONGER RUN
######################################################

#correlation_mat_file = wasp_visualization_dir + 'eqtl_correlation_mat_' + str(version) + '.txt'
#print_qtls(qtls, correlation_mat_file)


#raw_data = np.loadtxt(correlation_mat_file,dtype=str,delimiter='\t')
#test_ids = raw_data[:,0]
#pvalues = raw_data[:,1:].astype(float)

#correlation_pearson_mat = np.corrcoef(np.transpose(pvalues))
#correlation_spearman_mat = scipy.stats.spearmanr(pvalues)[0]


#fig = plt.figure()
#sns.heatmap(correlation_pearson_mat,cbar_kws={'label': 'Pearson Correlation'})
#plt.xlabel('Time Step')
#plt.ylabel('Time Step')
#output_file = wasp_visualization_dir + 'eqtl_pearson_correlation_heatmap_' + str(version) + '.png'
#fig.savefig(output_file)


#fig = plt.figure()
#sns.heatmap(correlation_spearman_mat,cbar_kws={'label': 'Spearman Correlation'})
#plt.xlabel('Time Step')
#plt.ylabel('Time Step')
#output_file = wasp_visualization_dir + 'eqtl_spearman_correlation_heatmap_' + str(version) + '.png'
#fig.savefig(output_file)


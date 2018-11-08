import numpy as np 
import os
import sys
import scipy.stats
import pdb
import gzip





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


def extract_variant_gene_pairs(cht_output_dir, parameter_string,efdr_thresh):
    qtls = {}
    for time_step in range(16):
        file_name = cht_output_dir + 'cht_results_' + parameter_string + '_time_' + str(time_step) + '_efdr_thresh_' + efdr_thresh + '_significant.txt'
        qtls = extract_all_significant_variant_gene_pairs(qtls, file_name)

    qtls = fill_in_pvalues(cht_output_dir, parameter_string, qtls)
    return qtls

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








def extract_egenes_with_smallest_geometric_mean_across_time(cht_output_dir, parameter_string,efdr_thresh):
    # Extract dictionary called qtls with:
    ## 1. keys that are all $gene_id + '_' + $variant_id (test name) pairs that are genome wide significant in 1 time step
    ## 2. values are an array of length number of time steps. Where each element in the arrray is the pvalue of that test in that time step
    qtls = extract_variant_gene_pairs(cht_output_dir, parameter_string, efdr_thresh)
    # Create dictionary called egenes with keys:
    ## 1. that are all $gene_id + '_' + $variant_id (test name) pairs such that the variant has the smallest geometric mean pvalue of all other variants
    egenes = get_egenes_through_geometric_mean(qtls)
    return egenes

def extract_gtex_file_names(used_gtex_tissues_file):
    tissues = []
    eqtl_files = []
    head_count = 0
    f = open(used_gtex_tissues_file)
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        tissues.append(data[0])
        eqtl_files.append(data[2])
    return np.asarray(tissues), np.asarray(eqtl_files)

def fill_in_per_time_step_eqtls(per_time_step_eqtl_pvalues, variant_gene_pair_mapping, cht_output_dir, parameter_string):
    valid_rows = {}
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
            pvalue = float(data[8])
            if test_id in variant_gene_pair_mapping:
                line_num = variant_gene_pair_mapping[test_id]
                per_time_step_eqtl_pvalues[line_num, time_step] = pvalue
                if line_num not in valid_rows:
                    valid_rows[line_num] = 1
                else:
                    valid_rows[line_num] = valid_rows[line_num] + 1
        f.close()
    return per_time_step_eqtl_pvalues, valid_rows

def fill_in_gtex_eqtls(gtex_eqtl_pvalues, variant_gene_pair_mapping, snp_to_rsid, all_eqtl_test_files):
    valid_rows = {}
    for i, tissue_file in enumerate(all_eqtl_test_files):
        print(i)
        f = gzip.open(tissue_file)
        head_count = 0
        for line in f:
            line = line.rstrip()
            data = line.split()
            if head_count == 0:
                head_count = head_count + 1
                continue
            snp_id = data[1]
            if snp_id not in snp_to_rsid:
                continue
            rs_id = snp_to_rsid[snp_id]
            ensamble_id = data[0].split('.')[0]
            test_id = ensamble_id + '_' + rs_id
            pvalue = float(data[6])
            if test_id in variant_gene_pair_mapping:
                line_num = variant_gene_pair_mapping[test_id]
                gtex_eqtl_pvalues[line_num, i] = pvalue
                if line_num not in valid_rows:
                    valid_rows[line_num] = 1
                else:
                    valid_rows[line_num] = valid_rows[line_num] + 1
        f.close()
    return gtex_eqtl_pvalues, valid_rows


def fill_in_ipsc_cm_eqtls(ipsc_cm_eqtl_pvalues, variant_gene_pair_mapping, cm_eqtl_file):
    valid_rows = {}
    head_count = 0
    f = open(cm_eqtl_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[0]
        rs_id = data[1]
        pvalue = float(data[2])
        test_id = ensamble_id + '_' + rs_id
        if test_id in variant_gene_pair_mapping:
            line_num = variant_gene_pair_mapping[test_id]
            ipsc_cm_eqtl_pvalues[line_num, 0] = pvalue
            valid_rows[line_num] = 1
    f.close()
    return ipsc_cm_eqtl_pvalues, valid_rows

def fill_in_ipsc_eqtls(ipsc_eqtl_pvalues, variant_gene_pair_mapping, ipsc_eqtl_file):
    valid_rows = {}
    f = open(ipsc_eqtl_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        ensamble_id = data[0].split('.')[0]
        rs_id = data[1].split('.')[0]
        test_id = ensamble_id + '_' + rs_id
        pvalue = float(data[-2])
        if test_id in variant_gene_pair_mapping:
            line_num = variant_gene_pair_mapping[test_id]
            ipsc_eqtl_pvalues[line_num, 0] = pvalue
            valid_rows[line_num] = 1
    f.close()
    return ipsc_eqtl_pvalues, valid_rows


def merge_time_steps_eqtls_and_banovich_eqtls(output_file, gene_variant_pairs, parameter_string, cht_output_dir, cm_eqtl_file, ipsc_eqtl_file):
    ordered_variant_gene_pairs = sorted(gene_variant_pairs.keys())
    variant_gene_pair_mapping = {}
    for i,val in enumerate(ordered_variant_gene_pairs):
        variant_gene_pair_mapping[val] = i
    
    num_rows = len(gene_variant_pairs)

    # IPSC eqtl matrix
    ipsc_eqtl_pvalues = np.zeros((num_rows, 1))
    ipsc_eqtl_pvalues, ipsc_eqtl_valid_rows = fill_in_ipsc_eqtls(ipsc_eqtl_pvalues, variant_gene_pair_mapping, ipsc_eqtl_file)

    # CM eqtl matrix
    ipsc_cm_eqtl_pvalues = np.zeros((num_rows,1))
    ipsc_cm_eqtl_pvalues, ipsc_cm_valid_rows = fill_in_ipsc_cm_eqtls(ipsc_cm_eqtl_pvalues, variant_gene_pair_mapping, cm_eqtl_file)

    #Per-time-step matrix
    per_time_step_eqtl_pvalues = np.zeros((num_rows, 16))
    per_time_step_eqtl_pvalues, per_time_step_valid_rows = fill_in_per_time_step_eqtls(per_time_step_eqtl_pvalues, variant_gene_pair_mapping, cht_output_dir, parameter_string)
    
    t = open(output_file, 'w')
    t.write('test_id')
    for time_step in range(16):
        t.write('\t' + 'time_step_' + str(time_step))
    t.write('\t' + 'ipsc')
    t.write('\t' + 'ipsc_cm')
    t.write('\n')

    for row_number, test_name in enumerate(ordered_variant_gene_pairs):
        if row_number in per_time_step_valid_rows and row_number in ipsc_cm_valid_rows and row_number in ipsc_eqtl_valid_rows:
            if per_time_step_valid_rows[row_number] == 16 and ipsc_cm_valid_rows[row_number] == 1 and ipsc_eqtl_valid_rows[row_number] == 1:
                t.write(test_name + '\t')
                t.write('\t'.join(per_time_step_eqtl_pvalues[row_number,:].astype(str)) + '\t')
                t.write(str(ipsc_eqtl_pvalues[row_number,0]) + '\t')
                t.write(str(ipsc_cm_eqtl_pvalues[row_number,0]) + '\n')
    t.close()



def merge_time_steps_eqtls_and_gtex_eqtls(output_file, gene_variant_pairs, parameter_string, cht_output_dir, used_gtex_tissues_file, snp_id_to_rs_id_file, snp_to_rsid, cm_eqtl_file, ipsc_eqtl_file):
    tissues, all_eqtl_test_files = extract_gtex_file_names(used_gtex_tissues_file)

    ordered_variant_gene_pairs = sorted(gene_variant_pairs.keys())
    variant_gene_pair_mapping = {}
    for i,val in enumerate(ordered_variant_gene_pairs):
        variant_gene_pair_mapping[val] = i
    
    num_rows = len(gene_variant_pairs)

    # IPSC eqtl matrix
    ipsc_eqtl_pvalues = np.zeros((num_rows, 1))
    ipsc_eqtl_pvalues, ipsc_eqtl_valid_rows = fill_in_ipsc_eqtls(ipsc_eqtl_pvalues, variant_gene_pair_mapping, ipsc_eqtl_file)

    # CM eqtl matrix
    ipsc_cm_eqtl_pvalues = np.zeros((num_rows,1))
    ipsc_cm_eqtl_pvalues, ipsc_cm_valid_rows = fill_in_ipsc_cm_eqtls(ipsc_cm_eqtl_pvalues, variant_gene_pair_mapping, cm_eqtl_file)

    # GTEx eqtl matrix
    gtex_eqtl_pvalues = np.zeros((num_rows,len(tissues)))
    gtex_eqtl_pvalues, gtex_valid_rows = fill_in_gtex_eqtls(gtex_eqtl_pvalues, variant_gene_pair_mapping, snp_to_rsid, all_eqtl_test_files)

    #Per-time-step matrix
    per_time_step_eqtl_pvalues = np.zeros((num_rows, 16))
    per_time_step_eqtl_pvalues, per_time_step_valid_rows = fill_in_per_time_step_eqtls(per_time_step_eqtl_pvalues, variant_gene_pair_mapping, cht_output_dir, parameter_string)
    
    t = open(output_file, 'w')
    t.write('test_id')
    for time_step in range(16):
        t.write('\t' + 'time_step_' + str(time_step))
    for tissue in tissues:
        t.write('\t' + tissue)
    t.write('\t' + 'ipsc')
    t.write('\t' + 'ipsc_cm')
    t.write('\n')

    for row_number, test_name in enumerate(ordered_variant_gene_pairs):
        if row_number in gtex_valid_rows and row_number in per_time_step_valid_rows and row_number in ipsc_cm_valid_rows and row_number in ipsc_eqtl_valid_rows:
            if gtex_valid_rows[row_number] == len(tissues) and per_time_step_valid_rows[row_number] == 16 and ipsc_cm_valid_rows[row_number] == 1 and ipsc_eqtl_valid_rows[row_number] == 1:
                t.write(test_name + '\t')
                t.write('\t'.join(per_time_step_eqtl_pvalues[row_number,:].astype(str)) + '\t')
                t.write('\t'.join(gtex_eqtl_pvalues[row_number,:].astype(str)) + '\t')
                t.write(str(ipsc_eqtl_pvalues[row_number,0]) + '\t')
                t.write(str(ipsc_cm_eqtl_pvalues[row_number,0]) + '\n')
    t.close()


def get_snp_id_to_rs_id_mapping(snp_id_to_rs_id_file, dosage_genotype_file):
    rs_ids = {}
    f = open(dosage_genotype_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#'):
            continue
        rs_id = data[2]
        if rs_id == '.':
            continue
        chrom = data[0]
        pos = data[1]
        ref = data[3]
        alt = data[4]
        rs_ids[rs_id] = (chrom,pos,ref,alt)
    f.close()
    snp_to_rs = {}
    county = 0
    f = open(snp_id_to_rs_id_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        rs_id = data[6]
        snp_id = data[2]
        ref = data[3]
        alt = data[4]
        chrom = data[0] 
        pos = data[1]
        if rs_id in rs_ids:
            tupler = rs_ids[rs_id]
            # Make sure its the same variant-gene pair
            if tupler[0] == data[0] and tupler[1] == data[1] and tupler[2] == data[3] and tupler[3] == data[4]:
                snp_to_rs[snp_id] = rs_id
            else:
                county=county + 1
    return snp_to_rs

def extract_all_egenes(cht_output_dir, parameter_string, cm_egene_file, ipsc_egene_file, used_gtex_tissues_file):
    egenes = {}
    #####################
    # Add ipsc data
    #####################
    f = open(ipsc_egene_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[0].split('.')[0]
        rs_id = data[1].split('.')[0]
        test_name = ensamble_id + '_' + rs_id
        egenes[test_name] = 1
    f.close()
    #####################
    # Add ipsc-cm data
    #####################
    f = open(cm_egene_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[0]
        rs_id = data[1]
        test_name = ensamble_id + '_' + rs_id
        egenes[test_name] = 1
    f.close()

    #####################
    # Add ipsc-temporal data
    #####################
    for time_step in range(16):
        time_step_egene_file = cht_output_dir + 'cht_results_' + parameter_string + '_time_' + str(time_step) + '_efdr_thresh_.05_significant_egenes.txt'
        f = open(time_step_egene_file)
        head_count = 0
        for line in f:
            line = line.rstrip()
            data = line.split()
            if head_count == 0:
                head_count = head_count + 1
                continue
            ensamble_id = data[1]
            rs_id = data[3]
            test_name = ensamble_id + '_' + rs_id
            egenes[test_name] = 1
        f.close()

    #####################
    # Add gtex egenes
    #####################
    f = open(used_gtex_tissues_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        egene_file = data[1]
        g = open(egene_file)
        head_count2 = 0
        for line2 in g:
            data2 = line2.rstrip().split()
            if head_count2 == 0:
                head_count2 = head_count2 + 1
                continue
            ensamble_id = data2[0].split('.')[0]
            rs_id = data2[18]
            test_name = ensamble_id + '_' + rs_id
            qvalue = float(data2[-2])
            if qvalue <= .01:
                egenes[test_name] = 1
        g.close()
    f.close()
    return egenes

def summarize_time_step_banovich_comparison(output_file, num_pc):
    aa = np.loadtxt(output_file,dtype=str,delimiter='\t')
    data = aa[1:,1:].astype(float)
    ipsc = data[:,-2]
    ipsc_cm = data[:,-1]
    ipsc_corrs = []
    ipsc_cm_corrs = []
    ipsc_log10_pvalues = []
    ipsc_cm_log10_pvalues = []
    for time_step in range(16):
        per_time_step_vec = data[:,time_step]
        ipsc_corrs.append(scipy.stats.spearmanr(ipsc, per_time_step_vec)[0])
        ipsc_cm_corrs.append(scipy.stats.spearmanr(ipsc_cm, per_time_step_vec)[0])

        ipsc_log10_pvalues.append(-np.log10(scipy.stats.spearmanr(ipsc, per_time_step_vec)[1]))
        ipsc_cm_log10_pvalues.append(-np.log10(scipy.stats.spearmanr(ipsc_cm, per_time_step_vec)[1]))
    #print(ipsc_corrs)
    time = range(16)
    plt.plot(time, ipsc_corrs, 'g')
    plt.plot(time, ipsc_cm_corrs, 'r')
    plt.savefig('joint_ipsc_ipsc_cm_correlation_' + str(num_pc) + '.png')

def extract_top_n_genes_from_a_time_step_eqtl_file(all_eqtl_file, num_genes):
    f = open(all_eqtl_file)
    head_count = 0
    gene_dicti = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[1]
        rs_id = data[3]
        pvalue = float(data[-1])
        if ensamble_id not in gene_dicti:
            gene_dicti[ensamble_id] = (rs_id, pvalue)
        else:
            old_pvalue = gene_dicti[ensamble_id][1]
            if pvalue < old_pvalue:
                gene_dicti[ensamble_id] = (rs_id, pvalue)
    f.close()
    arr = []
    for ensamble_id in gene_dicti.keys():
        tupler = gene_dicti[ensamble_id]
        rs_id = tupler[0]
        pvalue = tupler[1]
        arr.append((rs_id + '_' + ensamble_id, pvalue))
    arr.sort(key=lambda x: x[1])
    dicti = {}
    for i,val in enumerate(arr):
        if i < num_genes:
            dicti[val[0]] = 1
    return dicti

def find_pairs_in_ipsc_data(ipsc_eqtl_file, dicti):
    new_dicti = {}
    pvalue_arr = []
    f = open(ipsc_eqtl_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        ensamble_id = data[0].split('.')[0]
        rs_id = data[1].split('.')[0]
        if rs_id + '_' + ensamble_id in dicti:
            pvalue = float(data[-2])
            new_dicti[rs_id + '_' + ensamble_id] = pvalue 
            pvalue_arr.append(-np.log10(pvalue))
    f.close()
    return new_dicti

def find_pairs_in_ipsc_cm_data(cm_eqtl_file, dicti):
    new_dicti = {}
    pvalue_arr = []
    f = open(cm_eqtl_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[0]
        rs_id = data[1]
        if rs_id + '_' + ensamble_id in dicti:
            pvalue = float(data[2])
            new_dicti[rs_id + '_' + ensamble_id] = pvalue 
            pvalue_arr.append(-np.log10(pvalue + .000000000000001))
    f.close()
    return new_dicti

def extract_top_nn_from_each_time_step_and_banovich(num_genes, parameter_string, num_pc, cht_output_dir, cm_eqtl_file, ipsc_eqtl_file):
    output_file = cht_output_dir + parameter_string + '_top_' + str(num_genes) + '_eqtls_in_time_steps_for_banovich_results.txt'
    t = open(output_file, 'w')
    t.write('ensamble_id\trs_id\ttime_step\tdata_type\tpvalue\n')
    for time_step in range(16):
        print(time_step)
        all_eqtl_file = cht_output_dir + 'cht_results_' + parameter_string + '_time_' + str(time_step) + '_eqtl_results.txt'
        time_step_variant_gene_pair_dicti = extract_top_n_genes_from_a_time_step_eqtl_file(all_eqtl_file, num_genes)
        ipsc_eqtl_dicti = find_pairs_in_ipsc_data(ipsc_eqtl_file, time_step_variant_gene_pair_dicti)
        ipsc_cm_eqtl_dicti = find_pairs_in_ipsc_cm_data(cm_eqtl_file, time_step_variant_gene_pair_dicti)
        # PRINT
        for test_name in time_step_variant_gene_pair_dicti.keys():
            rs_id = test_name.split('_')[0]
            ensamble_id = test_name.split('_')[1]
            if test_name in ipsc_eqtl_dicti:
                ipsc_pvalue = ipsc_eqtl_dicti[test_name]
                t.write(ensamble_id + '\t' + rs_id + '\t' + str(time_step) + '\t' + 'iPSC' + '\t' + str(ipsc_pvalue) + '\n')
            if test_name in ipsc_cm_eqtl_dicti:
                ipsc_cm_pvalue = ipsc_cm_eqtl_dicti[test_name]
                t.write(ensamble_id + '\t' + rs_id + '\t' + str(time_step) + '\t' + 'iPSC_CM' + '\t' + str(ipsc_cm_pvalue) + '\n')
    t.close()


parameter_string = sys.argv[1]
cht_output_dir = sys.argv[2]
pc_num = sys.argv[3]
cm_eqtl_file = sys.argv[4]
ipsc_eqtl_file = sys.argv[5]
fdr = sys.argv[6]


parameter_string = parameter_string + '_num_pc_' + pc_num


####################################################
# Compare per time step eqtls with Nick's QTLs
####################################################
output_file = cht_output_dir + parameter_string + '_eqtls_across_time_steps_and_banovich_studies_geometric_mean_05.txt'
time_step_geometric_egenes_05 = extract_egenes_with_smallest_geometric_mean_across_time(cht_output_dir, parameter_string, fdr)
# Using these variant gene pairs, find those variant gene pairs in ipsc eqtl data and ipsc-cm eqtl data
merge_time_steps_eqtls_and_banovich_eqtls(output_file, time_step_geometric_egenes_05, parameter_string, cht_output_dir, cm_eqtl_file, ipsc_eqtl_file)





































#########################
# Old scripts (retired)
#########################



#snp_to_rsid = get_snp_id_to_rs_id_mapping(snp_id_to_rs_id_file, dosage_genotype_file)



####################################################
#  Find top nn genes (and variant) in each time step. 
# Then for each time step, get pvalues of those variant-gene pairs using nicks data
####################################################
#num_genes = 100
#extract_top_nn_from_each_time_step_and_banovich(num_genes, parameter_string, pc_num, cht_output_dir, cm_eqtl_file, ipsc_eqtl_file)

#num_genes = 200
#extract_top_nn_from_each_time_step_and_banovich(num_genes, parameter_string, pc_num, cht_output_dir, cm_eqtl_file, ipsc_eqtl_file)

#num_genes = 300
#extract_top_nn_from_each_time_step_and_banovich(num_genes, parameter_string, pc_num, cht_output_dir, cm_eqtl_file, ipsc_eqtl_file)



###########################################################################
# Method 1: Extract egenes with smallest geometric mean across time steps (efdr <= .05)
###########################################################################
#output_file = cht_output_dir + parameter_string + '_eqtls_across_time_steps_and_gtex_v7_geometric_mean_05.txt'
#time_step_geometric_egenes_05 = extract_egenes_with_smallest_geometric_mean_across_time(cht_output_dir, parameter_string, '.05')
# merge_time_steps_eqtls_and_gtex_eqtls(output_file, time_step_geometric_egenes_05, parameter_string, cht_output_dir, used_gtex_tissues_file, snp_id_to_rs_id_file, snp_to_rsid, cm_eqtl_file, ipsc_eqtl_file)


###########################################################################
# Method 2: Extract egenes with smallest geometric mean across time steps (efdr <= .1)
###########################################################################
#output_file = cht_output_dir + parameter_string + '_eqtls_across_time_steps_and_gtex_v7_geometric_mean_1.txt'
#time_step_geometric_egenes_1 = extract_egenes_with_smallest_geometric_mean_across_time(cht_output_dir, parameter_string, '.1')
#merge_time_steps_eqtls_and_gtex_eqtls(output_file, time_step_geometric_egenes_1, parameter_string, cht_output_dir, used_gtex_tissues_file, snp_id_to_rs_id_file, snp_to_rsid, cm_eqtl_file, ipsc_eqtl_file)


###########################################################################
# Method 3 Extract anything that is an egene in any gtex tissue, ipsc, ipsc-cm, or ipsc-time series data sets
###########################################################################
#output_file = cht_output_dir + parameter_string + '_eqtls_across_time_steps_and_gtex_v7_all_egenes.txt'
#all_egenes = extract_all_egenes(cht_output_dir, parameter_string, cm_egene_file, ipsc_egene_file, used_gtex_tissues_file)
#merge_time_steps_eqtls_and_gtex_eqtls(output_file, all_egenes, parameter_string, cht_output_dir, used_gtex_tissues_file, snp_id_to_rs_id_file, snp_to_rsid, cm_eqtl_file, ipsc_eqtl_file)



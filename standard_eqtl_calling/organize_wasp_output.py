import numpy as np
import os
import sys
import pdb
import gzip

# Create dictionary of all genes we are going to be testing (ie those that passed our filters)
def get_measured_genes(expression_mat):
    dicti = {}  # initialize
    head_count = 0  # For header
    f = open(expression_mat)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # skip header
            head_count = head_count +1
            continue
        # Gene ids are found in the first column of each line
        gene_id = data[0]
        dicti[gene_id] = -1
    return dicti

# Create dictionary mapping all genes found in measured_genes, and are located on chromosome $chrom_num, to their tss
def get_mapping_from_gene_to_tss(gencode_gene_annotation_file, chrom_num, measured_genes):
    gene_to_tss_mapping = {}
    f = gzip.open(gencode_gene_annotation_file)
    count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#'):  # ignore header lines
            continue
        gene_type = data[13].split('"')[1]  # ie protein_coding,pseudo_gene,etc
        gene_name = data[9].split('"')[1].split('.')[0]  # ensamble id
        line_chrom_num = data[0]
        start = int(data[3])  # Start  of gene
        end = int(data[4])  # End (downstream) of gene
        gene_part = data[2]  # gene,UTR,exon,etc
        strand = data[6]  # either positive or negative
        if line_chrom_num != 'chr' + chrom_num:  # limit to chromosome of interest
            continue
        if gene_name not in measured_genes:  # Only care about genes that we have measurements for
            continue
        # We now are limited to measured genes on the correct chromosome
        if gene_name not in gene_to_tss_mapping:  # We haven't seen this gene before
            if strand == '+':  # positive strand
                gene_to_tss_mapping[gene_name] = (start, strand)
            elif strand == '-':  # negative strand
                gene_to_tss_mapping[gene_name] = (end, strand)
        else:  # We've seen this gene before
            old_tuple = gene_to_tss_mapping[gene_name]
            old_tss = old_tuple[0]
            old_strand = old_tuple[1]
            tss = old_tss
            if old_strand != strand:  # Error checking
                print('ASSUMPTION ERROR')
                pdb.set_trace()
            if strand == '+' and start < old_tss: 
                tss = start
            elif strand == '-' and end > old_tss:
                tss = end
            gene_to_tss_mapping[gene_name] = (tss, strand)
    ordered_genes = []
    ordered_tss = []
    for gene_id in gene_to_tss_mapping.keys():
        tupler = gene_to_tss_mapping[gene_id]
        tss = tupler[0]
        ordered_genes.append(gene_id)
        ordered_tss.append(tss)
    return gene_to_tss_mapping, ordered_genes, ordered_tss

# Create dictionary mapping from gene target regions to gene names
def get_mapping_from_target_region_to_gene_name(target_region_file, chrom_num):
    f = open(target_region_file)
    dicti = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if data[0] != 'chr' + chrom_num:
            continue
        gene_id = data[6].split('_')[2]
        regions = data[7] + '_' + data[8]
        dicti[regions] = gene_id
    return dicti


def get_cell_lines_in_this_time_step(corrected_quantile_normalized_expression,time_step):
    head_count = 0
    dicti = {}
    f = open(corrected_quantile_normalized_expression)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            for ele in data[1:]:
                curr_time_step = ele.split('_')[1]
                curr_cell_line = ele.split('_')[0]
                if curr_time_step == time_step:
                    dicti[curr_cell_line] = 1
    return dicti

def get_maf(genotype_array):
    af = np.sum(genotype_array)/(2.0*len(genotype_array))
    return min(af, 1.0 - af)

def get_rsid_info(dosage_genotype_file, chrom_num, cell_lines):
    rsid_to_position = {}
    rsid_to_maf = {}
    f = open(dosage_genotype_file)
    indices = []
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#CHRO'):
            for i, cell_line in enumerate(data):
                if cell_line in cell_lines:
                    indices.append(i)
            indices = np.asarray(indices)
            continue
        if line.startswith('#'):
            continue
        if chrom_num != data[0]:
            continue
        pos = data[1]
        rsid = data[2]
        genotype = np.asarray(data)[indices]
        maf = get_maf(genotype.astype(float))
        rsid_to_position[rsid] = pos
        rsid_to_maf[rsid] = str(maf)
    return rsid_to_position, rsid_to_maf


def organize_results_driver(t, chrom_input_file, dosage_genotype_file, corrected_quantile_normalized_expression, gencode_gene_annotation_file, target_region_file, chrom_num, cell_lines):
    
    rsid_to_position, rsid_to_maf = get_rsid_info(dosage_genotype_file, str(chrom_num), cell_lines)

    # Create dictionary of all genes we are going to be testing (ie those that passed our filters)
    measured_genes = get_measured_genes(corrected_quantile_normalized_expression)
    # Create dictionary mapping all genes found in measured_genes, and are located on chromosome $chrom_num, to their tss
    gene_to_tss, ordered_genes, ordered_tss = get_mapping_from_gene_to_tss(gencode_gene_annotation_file, str(chrom_num), measured_genes)
    ordered_tss = np.asarray(ordered_tss)

    # Create dictionary mapping from gene target regions to gene names
    target_region_to_gene_names = get_mapping_from_target_region_to_gene_name(target_region_file, str(chrom_num))


    # Print to new file
    head_count = 0
    f = open(chrom_input_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        if data[0] != 'chr' + str(chrom_num):
            print('ERORORO!')
        rsid = data[2]
        regions = data[5] + '_' + data[6]

        rs_position = rsid_to_position[rsid]
        maf = rsid_to_maf[rsid]

        gene_id = target_region_to_gene_names[regions]

        gene_tss = str(gene_to_tss[gene_id][0])

        pvalue = data[10]
        alpha = data[11]
        beta = data[12]

        t.write(str(chrom_num) + '\t' + gene_id + '\t' + gene_tss + '\t' + rsid + '\t' + rs_position + '\t' + maf + '\t' + alpha + '\t' + beta + '\t' + pvalue + '\n')

    return t









time_step = sys.argv[1]
wasp_test_dir = sys.argv[2]
dosage_genotype_file = sys.argv[3]
corrected_quantile_normalized_expression = sys.argv[4]
gencode_gene_annotation_file = sys.argv[5]
target_region_file = sys.argv[6]
wasp_results_stem = sys.argv[7]


# Initialize output file
output_file = wasp_results_stem + 'eqtl_results.txt'
t = open(output_file, 'w')
t.write('chrom_num\tgene_id\tgene_position\tregulatory_site_id\tregulatory_site_position\tregulatory_allele_frequency\talpha\tbeta\tpvalue\n')

# Get list of cell lines that are used at this time step
cell_lines = get_cell_lines_in_this_time_step(corrected_quantile_normalized_expression,time_step)

# Loop through chromosomes
for chrom_num in range(1,23):
    # Chromosome specific input file
    chrom_input_file = wasp_results_stem + str(chrom_num) + '.txt'
    # Add this chromosome specific input file to our output file
    t = organize_results_driver(t, chrom_input_file, dosage_genotype_file, corrected_quantile_normalized_expression, gencode_gene_annotation_file, target_region_file, chrom_num, cell_lines)
t.close()
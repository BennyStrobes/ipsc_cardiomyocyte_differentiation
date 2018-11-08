import numpy as np
import os
import sys
import pdb
import gzip
import chromosome

import tables


########################################
# Code borrowed from bryce
#######################################
SNP_UNDEF = -1


class SNP(object):
    def __init__(self, chrom, pos, name, ref_allele, alt_allele):
        self.chrom = chrom
        self.pos = pos
        self.name = name
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        

class SNPFiles(object):
    def __init__(self, snp_tab,snp_index,haplotype):
        # open tracks where SNP information can be extracted
        self.snp_tab_h5 = tables.open_file(snp_tab, "r")
        self.snp_index_h5 = tables.open_file(snp_index, "r")
        self.hap_h5 = tables.open_file(haplotype, "r")


    def close(self):
        self.snp_tab_h5.close()
        self.snp_index_h5.close()
        self.hap_h5.close()
        


class CombinedFiles(object):
    def __init__(self, output_dir, chrom_list,time_step):
        time_step_str = str(time_step)
        # combined allele-specific read counts
        as_count_filename = "%s/combined_as_count_%s.h5" % (output_dir,time_step_str)
        self.as_count_h5 = tables.open_file(as_count_filename, "w")
        
        # combined mapped read counts
        read_count_filename = "%s/combined_read_count_%s.h5" % (output_dir,time_step_str)
        self.read_count_h5 = tables.open_file(read_count_filename, "w")

        # counts of genotypes
        ref_count_filename = "%s/combined_ref_count_%s.h5" % (output_dir,time_step_str)
        self.ref_count_h5 = tables.open_file(ref_count_filename, "w")
        
        alt_count_filename = "%s/combined_alt_count_%s.h5" % (output_dir,time_step_str)
        self.alt_count_h5 = tables.open_file(alt_count_filename, "w")
        
        het_count_filename = "%s/combined_het_count_%s.h5" % (output_dir,time_step_str)
        self.het_count_h5 = tables.open_file(het_count_filename, "w")
        
        self.filenames = [as_count_filename, read_count_filename,
                          ref_count_filename, alt_count_filename,
                          het_count_filename]

        self.h5_files = [self.as_count_h5, self.read_count_h5,
                         self.ref_count_h5, self.alt_count_h5, 
                         self.het_count_h5]

        # initialize all of these files
        atom = tables.UInt16Atom(dflt=0)
        
        for h5f in self.h5_files:
            for chrom in chrom_list:
                self.create_carray(h5f, chrom, atom)
                
    def create_carray(self, h5f, chrom, atom):
        zlib_filter = tables.Filters(complevel=1, complib="zlib")

        # create CArray for this chromosome
        shape = [chrom.length]
        carray = h5f.create_carray(h5f.root, chrom.name,
                                  atom, shape, filters=zlib_filter)

        return carray


    
    def add_counts(self, chrom_list, ind_files, snp_files, ind_idx):
        # add contribution from one individual to combined counts
        for chrom in chrom_list:
            node_name = "/%s" % chrom.name

            if node_name not in ind_files.ref_as_count_h5:
                continue
            if node_name not in ind_files.alt_as_count_h5:
                continue
            if node_name not in ind_files.read_count_h5:
                continue
            if node_name not in snp_files.hap_h5:
                continue
            if node_name not in snp_files.snp_index_h5:
                continue
            if node_name not in snp_files.snp_tab_h5:
                continue

            sys.stderr.write("  %s\n" % chrom.name)
            
            node = self.as_count_h5.get_node(node_name)
                
            ind_node_ref = ind_files.ref_as_count_h5.get_node(node_name)
            ind_node_alt = ind_files.alt_as_count_h5.get_node(node_name)
            
            node[:] += np.minimum(ind_node_alt[:], ind_node_ref[:])

            node = self.read_count_h5.get_node(node_name)
            ind_node = ind_files.read_count_h5.get_node("/%s" % chrom.name)
            node[:] += ind_node[:]


            # get haplotypes for this individual
            hap_a_idx = ind_idx*2
            hap_b_idx = ind_idx*2+1            
            hap_tab = snp_files.hap_h5.get_node("/%s" % chrom.name)
            a_hap = hap_tab[:, hap_a_idx]
            b_hap = hap_tab[:, hap_b_idx]

            # determine genotype of SNPs for this individual
            is_homo_ref = (a_hap == 0) & (b_hap == 0)
            is_het      = ((a_hap == 0) & (b_hap == 1)) | ((a_hap == 1) & (b_hap == 0))
            is_homo_alt = (a_hap == 1) & (b_hap == 1)


            # get genomic location of SNPs
            i = snp_files.snp_index_h5.get_node("/%s" % chrom.name)[:]
            chrom_idx = np.where(i != SNP_UNDEF)[0]
            snp_idx = i[chrom_idx]
  

            # add to total genotype counts
            node = self.ref_count_h5.get_node("/%s" % chrom.name)
            node[chrom_idx] += is_homo_ref[snp_idx]
            
            node = self.het_count_h5.get_node("/%s" % chrom.name)
            node[chrom_idx] += is_het[snp_idx]
            
            node = self.alt_count_h5.get_node("/%s" % chrom.name)
            node[chrom_idx] += is_homo_alt[snp_idx]
    
            

    def close(self):
        for h5f in self.h5_files:
            h5f.close()

        # remove tempory combined files
        for filename in self.filenames:
            os.unlink(filename)
                         



class CountFiles(object):    
    def __init__(self, read_count_dir, individual):
        # open read count tracks for a single individual
        self.ref_as_count_h5 = tables.open_file("%s/ref_as_counts.%s.h5" % 
                                            (read_count_dir, individual), "r")
        self.alt_as_count_h5 = tables.open_file("%s/alt_as_counts.%s.h5" % 
                                            (read_count_dir, individual), "r")
        self.read_count_h5 = tables.open_file("%s/read_counts.%s.h5" %
                                             (read_count_dir, individual), "r")

        
    def close(self):
        """closes all of the data files"""
        self.ref_as_count_h5.close()
        self.alt_as_count_h5.close()
        self.read_count_h5.close()
        



def get_samples_index(file_name):
    """Gets dictionary of sample_id => index mappings that is used 
    to lookup information in the genotype and haplotype tables"""
    f = open(file_name)

    ind_dict = {}
    
    idx = 0
    for line in f:
        words = line.rstrip().split()

        # name = words[0].replace("NA", "")
        name = words[0]

        if name in ind_dict:
            raise ValueError("sample identifier '%s' appears multiple "
                             "times in file %s" % (name, options.samples))
        
        ind_dict[name] = idx
                
        idx += 1

    return ind_dict

def get_individual_array(file_name):
    arr = []
    f = open(file_name)
    for line in f:
        line = line.rstrip()
        arr.append(line)
    return arr

def extract_rna_info(chrom_info_file, raw_allelic_counts_dir, genotype_dir,time_step, target_regions_dir):
    # make dictionary of identifier => index mapping
    all_genotype_samples_file = genotype_dir + 'all_genotyped_samples.txt'
    samp_idx = get_samples_index(all_genotype_samples_file)

    # Initialize chromosome objects
    chrom_list = chromosome.get_all_chromosomes(chrom_info_file)
    chrom_dict = chromosome.get_chromosome_dict(chrom_info_file)

    snp_files = SNPFiles(genotype_dir + 'snp_tab.h5',genotype_dir + 'snp_index.h5',genotype_dir+'haps.h5')

    # STEP 1: make combined HDF5 files of AS counts, 
    # total mapped read counts, and genotype counts
    individuals = get_individual_array(target_regions_dir + 'rna_seq_samples_' + str(time_step) + '.txt')
    combined_files = CombinedFiles(raw_allelic_counts_dir, chrom_list,time_step)
    for ind in individuals:
        print(ind)
        sample_id = ind + '_' + str(time_step)
        count_files = CountFiles(raw_allelic_counts_dir, sample_id)
            
        ind_idx = samp_idx[ind]
        combined_files.add_counts(chrom_list, count_files, snp_files, ind_idx)

        count_files.close()

    return combined_files





########################################
# my code
#######################################



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


def get_non_overlapping_exons(arr):
    chromosome = np.zeros(259250621)
    start = 259250650
    end = 0
    for exon in arr:
        exon_start = exon[0]
        exon_end = exon[1]
        if exon_start < start:
            start = exon_start
        if exon_end > end:
            end = exon_end
        chromosome[exon_start:(exon_end+1)] = np.ones(exon_end-exon_start+1)
    non_overlapping_array = []

    in_exon = False
    for pos in range(start-1,end+2):
        if in_exon == False and chromosome[pos] == 1.0:  # Now starting an exon
            temp_start = pos
            in_exon = True
        elif in_exon == True and chromosome[pos] == 0.0: # Just left an exon
            temp_end = pos - 1
            in_exon = False
            non_overlapping_array.append((temp_start,temp_end))
    if in_exon == True:
        print('assumption error')
        pdb.set_trace()
    return non_overlapping_array


def get_mapping_from_gene_to_exonic_region(gencode_gene_annotation_file, chrom_num, gene_to_tss):
    gene_to_exonic = {}
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
        if gene_name not in gene_to_tss:  # Only care about genes that we have measurements for
            continue
        if gene_part != 'exon':  # Limit to exonic regions
            continue
        if end < start:
            print('FATAL ASSUMPTION ERROR')
            pdb.set_trace()

        if gene_name not in gene_to_exonic:
            gene_to_exonic[gene_name] = []
        gene_to_exonic[gene_name].append((start,end))

    gene_to_exonic_non_overlapping = {}
    for i, gene_id in enumerate(gene_to_exonic.keys()):
        arr = gene_to_exonic[gene_id]
        non_overlapping_array = get_non_overlapping_exons(arr)
        gene_to_exonic_non_overlapping[gene_id] = non_overlapping_array
    return gene_to_exonic_non_overlapping

# Extract list of cell lines that are found in all time steps (some time steps have different numbers of cell lines)
# When calling maf cutoff, we are to use only cell lines found in all time steps
def extract_list_of_cell_lines_in_all_time_steps(quantile_normalized_expression, time_step):
    # First, extract array of sample ids from quantile normalized expression matrix
    f = open(quantile_normalized_expression)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            sample_ids = data[1:]
            continue
        continue
    f.close()
    # Now loop through all time steps to determine which cell lines were observed in which time step
    cell_lines_all_time_steps = []
    temp_dicti = {}
    for sample_id in sample_ids:
        # check if sample_id is from the current iteration's time step
        if sample_id.endswith('_' + str(time_step)) == False:
            continue
        # Get cell line from sample id
        cell_line = sample_id.split('_')[0]
        temp_dicti[cell_line] = 1
    return temp_dicti

# Extract an array of length number of time steps, where each element of the array is another array that contains indices of cell lines observed for this time step
def get_indices_for_each_time_step(data, cell_lines_all_time_steps):
    indices = []
    for i, val in enumerate(data):
        if val in cell_lines_all_time_steps:
            indices.append(i)
    return indices

def get_maf(genotype_array):
    af = np.sum(genotype_array)/(2.0*len(genotype_array))
    return min(af, 1.0 - af)


def pass_maf_cutoff_all_time_steps(data,indices_all_time_steps, maf_cutoff):
    pass_filter = True
    for time_step,index_vector in enumerate(indices_all_time_steps):
        genotype_for_one_time_step = data[index_vector].astype(float)
        if get_maf(genotype_for_one_time_step) < maf_cutoff:
            pass_filter = False
    return pass_filter

def region_to_string(region):
    starts = []
    ends = []
    for exon in region:
        starts.append(exon[0])
        ends.append(exon[1])
    starts = np.asarray(starts).astype(str)
    ends = np.asarray(ends).astype(str)
    return ';'.join(starts), ';'.join(ends)


def pass_read_count_filter(start_string, end_string, read_counts, as_read_counts, min_read_count,min_as_read_count, het_counts, min_het_count):
    starts = np.asarray(start_string.split(';')).astype(int)
    ends = np.asarray(end_string.split(';')).astype(int)

    total_reads = 0
    total_as_reads = 0
    total_hets = 0
    for index, start in enumerate(starts):
        end = ends[index]
        n_reads = np.sum(read_counts[start-1:end])
        n_as_reads = np.sum(as_read_counts[start-1:end])
        n_hets = np.sum(het_counts[start-1:end])
        total_reads = total_reads + n_reads
        total_as_reads = total_as_reads + n_as_reads
        total_hets = total_hets + n_hets
    if (total_reads >= min_read_count) and (total_as_reads >= min_as_read_count) and (total_hets >= min_het_count):
        bin = True
    else:
        bin = False
    return bin

def get_target_regions(chrom_num, gene_expression_input_file, gencode_gene_annotation_file, dosage_genotype_file, cis_distance, maf_cutoff,t, min_read_count, min_as_read_count, time_step, read_counts, as_read_counts, het_counts, min_het_count):
    # Create dictionary of all genes we are going to be testing (ie those that passed our filters)
    measured_genes = get_measured_genes(gene_expression_input_file)
    # Create dictionary mapping all genes found in measured_genes, and are located on chromosome $chrom_num, to their tss
    gene_to_tss, ordered_genes, ordered_tss = get_mapping_from_gene_to_tss(gencode_gene_annotation_file, str(chrom_num), measured_genes)
    ordered_tss = np.asarray(ordered_tss)
    # Create dictionary mapping all measured genes to list of their exonic regions
    gene_to_region = get_mapping_from_gene_to_exonic_region(gencode_gene_annotation_file, str(chrom_num), gene_to_tss)


    # Extract list of length number of time steps where each element is a dictionary that contains the cell lines observed for that time step
    # When calling maf cutoff, we require it to pass the maf in each observed time step
    cell_lines_all_time_steps = extract_list_of_cell_lines_in_all_time_steps(gene_expression_input_file, time_step)

    f = open(dosage_genotype_file)
    str_chrom_num = str(chrom_num)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#CHROM'):  # HEADER with sample names
            # Extract an array of length number of time steps, where each element of the array is another array that contains indices of cell lines observed for this time step
            indices = get_indices_for_each_time_step(data, cell_lines_all_time_steps)
            continue
        if line.startswith('#'):  # headers we do not care about for now
            continue
        # normal line
        # Extract relavent features of line
        rs_id = data[2]
        chromer = 'chr' + data[0]
        pos = data[1]
        ref_allele = data[3]
        alt_allele = data[4]
        if str_chrom_num != data[0]:  # Not on correct chromsome
            continue

        if rs_id == '.':  #remove variants that do not have an rsID
            continue

        # Ignore snps that have maf < cutoff in any time step
        genotype_array = np.asarray(data)[indices].astype(float)
        if get_maf(genotype_array) < maf_cutoff:
            continue

        distance_to_tss = abs(ordered_tss - int(pos))

        # Get indices of all genes that have a TSS within cis_distance of this variant
        genes_in_cis_distance = np.where(distance_to_tss <= cis_distance)[0]

        # Loop through indices of all these genes genes
        for index in genes_in_cis_distance:
            gene_id = ordered_genes[index]
            gene_tss = ordered_tss[index]
            region = gene_to_region[gene_id]
            start_string, end_string = region_to_string(region)

            if pass_read_count_filter(start_string, end_string, read_counts, as_read_counts, min_read_count,min_as_read_count, het_counts, min_het_count) == True:
                t.write(chromer + ' ' + str(pos) + ' ' + str(int(pos)+1) + ' ' + ref_allele + ' ' + alt_allele + ' + ' + chromer + '_' + rs_id + '_' + gene_id + ' ' + start_string + ' ' + end_string + '\n')
    return t



gene_expression_input_file = sys.argv[1]
gencode_gene_annotation_file = sys.argv[2]
dosage_genotype_file = sys.argv[3]
genotype_dir = sys.argv[4]
cis_distance = int(sys.argv[5])
maf_cutoff = float(sys.argv[6])
raw_allelic_counts_dir = sys.argv[7]
min_read_count = int(sys.argv[8])
min_as_read_count = int(sys.argv[9])
chrom_info_file = sys.argv[10]
time_step = sys.argv[11]
min_het_count = int(sys.argv[12])
target_regions_dir = sys.argv[13]

# Extract data structure that contains number of raw read counts at a specific genomic position, summed across all the samples
combined_files = extract_rna_info(chrom_info_file, raw_allelic_counts_dir, genotype_dir, time_step, target_regions_dir)

# Open output file handle
t = open(target_regions_dir + 'target_regions_cis_distance_' + str(cis_distance) + '_maf_cutoff_' + str(maf_cutoff) + '_min_reads_' + str(min_read_count) + '_min_as_reads_' + str(min_as_read_count) + '_min_het_counts_' + str(min_het_count) + '_time_' + str(time_step) +'.txt','w')

# Loop through chromosomes
for chrom_num in range(1,23):
    print(chrom_num)

    node_name = '/chr' + str(chrom_num)
    # Total number of read counts (summed across all samples) for all positions on this chromosome
    read_counts = combined_files.read_count_h5.get_node(node_name)[:]
    # Number of read counts (summed across all samples) on less popular allele (min transformation) for all positions on this chromosome
    as_read_counts = combined_files.as_count_h5.get_node(node_name)[:]
    # Number of heterozygous sites (summed across all samples) for all positions on this chromosome
    het_counts = combined_files.het_count_h5.get_node(node_name)[:]


    t = get_target_regions(chrom_num, gene_expression_input_file, gencode_gene_annotation_file, dosage_genotype_file, cis_distance, maf_cutoff,t, min_read_count, min_as_read_count, time_step, read_counts, as_read_counts, het_counts, min_het_count)


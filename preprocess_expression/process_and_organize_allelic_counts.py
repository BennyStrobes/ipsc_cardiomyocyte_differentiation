import numpy as np
import os
import sys
import pdb
import gzip


# Not all samples are to be used. We performed filtering on samples in the total expression preprocessing.
# Retrieves ordered list of samples to be used (based on total expression filtering)
def get_ordered_samples_from_total_expression(file_name):
    f = open(file_name)
    head_count = 0  # Indicator for header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # If header
            head_count = head_count + 1
            ordered_samples = np.asarray(data[1:])
        else:  # If we've gone past the header, break out of loop
            break
    return ordered_samples


# Create a dictionary whos keys are sites and values are another (nested dictionary). The keys of the nested dictionary are sample_ids, and the values of the nested dictionary are  read counts of the form refAlleleCounts_totalCounts
# So the dictionary can be accessed as follows dictionary[site_id][sample_id] to get the number of reads mapping for site_id in sample_id
def update_allelic_counts_dictionary(allelic_counts_dictionary, sample_specific_read_count_file, sample_id):
    f = open(sample_specific_read_count_file)
    head_count = 0  # Used to skip header
    #  Loop through input file (each line is a potential heterozygous site)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        # Create unique string identifier for each site (of the form chromNum_position_rsID_refAllele_altAllele)
        site_id = data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3] + '_' + data[4]
        #  Create string corresponding to number of read counts at (sample, site) pair.
        #  String is of the form refAlleleCounts_totalCounts
        counts_string = data[5] + '_' + data[7]

        #  Add to dictionary
        if site_id not in allelic_counts_dictionary:  # If site never seen before, add site to dictoinary
            allelic_counts_dictionary[site_id] = {}

        # ERROR CHECKING
        if sample_id in allelic_counts_dictionary[site_id]:
            print('Fatal error: Saw sampleID,site pair twice')
            pdb.set_trace()
        #  Add counts to dictionary
        allelic_counts_dictionary[site_id][sample_id] = counts_string
    return allelic_counts_dictionary


# Helper function to convert from dictionary to ordered array
def extract_ordered_array_of_counts(dictionary, ordered_samples):
    # Initialize count_array to be all Nan
    count_array = ['Nan']*len(ordered_samples)
    # Loop through all samples to check if that sample has reads mapped to this site
    for i, sample_id in enumerate(ordered_samples):
        # Check if sample has reads mapped
        if sample_id in dictionary:
            # If it does, update count_array
            count_array[i] = dictionary[sample_id]
    return np.asarray(count_array)


# Main driver for part 1: Merge all of the sample specific count files into one table
def merge_sample_specific_count_files(raw_allelic_counts_dir, ordered_samples, raw_allelic_counts_output_file):
    # Create a dictionary whos keys are sites and values are another (nested dictionary). The keys of the nested dictionary are sample_ids, and the values of the nested dictionary are read counts of the form refAlleleCounts_totalCounts
    # So the dictionary can be accessed as follows dictionary[site_id][sample_id] to get the number of reads mapping for site_id in sample_id
    allelic_counts_dictionary = {}
    for sample_id in ordered_samples:
        sample_specific_read_count_file = raw_allelic_counts_dir + sample_id + '_merged_wasp_corrected.txt'  # output from GATK ASEReadCounter
        allelic_counts_dictionary = update_allelic_counts_dictionary(allelic_counts_dictionary, sample_specific_read_count_file, sample_id)

    # Now we want to print the counts to one ouptut file
    t = open(raw_allelic_counts_output_file, 'w')
    #  Write header
    t.write('siteID\t' + '\t'.join(ordered_samples) + '\n')
    # Loop through all sites and print each site (line)
    for site_id in sorted(allelic_counts_dictionary.keys()):
        t.write(site_id + '\t')  # Row header
        count_array = extract_ordered_array_of_counts(allelic_counts_dictionary[site_id], ordered_samples)  # Helper function to convert from dictionary to ordered array
        t.write('\t'.join(count_array) + '\n')
    t.close()

#Fill in chromosome object from chromosome[start:end+1] with the addition of this gene name
def add_gene_to_chromosome_object(chromosome,start,end,gene_name):
    for pos in range(start,end+1):
        if chromosome[pos] == 'NULL': #No genes already mapped to this position
            chromosome[pos] = gene_name
        else: #at least one gene already mapped to this position
            chromosome[pos] = chromosome[pos] + '_' + gene_name
    return chromosome


#  Create an array of length(chromosome) [in bp]. Value of an element corresponds to a list of genes that overlap that basepair. This is used for efficient searching.
def make_chromosome_with_gene_names(chromosome, chrom_num, gencode_gene_annotation_file):
    # loop through gencode file
    f = gzip.open(gencode_gene_annotation_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#'):  # ignore header lines
            continue
        gene_type = data[13].split('"')[1]  # ie protein_coding,pseudo_gene,etc
        gene_name = data[9].split('"')[1]  # ensamble id
        line_chrom_num = data[0]
        gene_part = data[2]  # gene,UTR,exon,etc
        if gene_type != 'protein_coding':  # limit to only protein coding genes
            continue
        if line_chrom_num != 'chr' + chrom_num:  # limit to chromosome of interest
            continue
        if gene_part != 'exon' and gene_part != 'UTR':  # Only consider exons and UTRs
            continue
        start = int(data[3])
        end = int(data[4])
        #Fill in chromosome object from chromosome[start:end+1] with the addition of this gene name
        chromosome = add_gene_to_chromosome_object(chromosome, start, end, gene_name)
    return chromosome

# Main driver for part 2: Map the sites to protein coding genes and remove sites that don't lie on a protein-coding gene
def map_sites_to_genes(raw_allelic_counts_output_file, allelic_counts_gene_mapped_output_file, gencode_gene_annotation_file):
     # open output  file file-handle
    t = open(allelic_counts_gene_mapped_output_file, 'w')
    # Write header to output file
    f = open(raw_allelic_counts_output_file)
    header_line = f.next().rstrip()
    t.write(header_line + '\n')
    f.close()  #  close input file handle

    # Map loci to genes for each chromosome seperately. So loop through chromosomes
    for chrom_num in range(1,23):
        print(chrom_num)
        #  Create an array of length(chromosome) [in bp]. Value of an element corresponds to a list of genes that overlap that basepair. This is used for efficient searching.
        chromosome = ['NULL']*259250621  # Initialize the array
        #  Fill in array with genes that overlap that position
        chromosome = make_chromosome_with_gene_names(chromosome, str(chrom_num), gencode_gene_annotation_file)

        # Stream previous allelic counts file, and print the line if it belongs to this chromosome and it overlaps a gene
        f = open(raw_allelic_counts_output_file)
        head_count = 0  # to be used to skip header
        chrom_stringer = str(chrom_num) + '_'  # To be used to check which chromosome site is on
        for line in f:
            line = line.rstrip()
            if head_count == 0:  # Skip header
                head_count = head_count + 1
                continue
            #  Check if site falls on desired chromsome
            if line.startswith(chrom_stringer) == False:  # Site is not on this chromosome
                continue
            #  Parse line a bit
            data = line.split()
            site_id = data[0]
            site_id_info = site_id.split('_')
            site_position = int(site_id_info[1])

            gene_stringer = '_'.join(np.unique(np.asarray(chromosome[site_position].split('_'))))  # Take unique genes
            if gene_stringer == 'NULL':  #  No genes overlapping site
                continue

            new_site_id = site_id + '_' + gene_stringer  # Add gene name to end of site_id
            t.write(new_site_id + '\t' + '\t'.join(data[1:]) + '\n')
        f.close()
    t.close()

#  Create list of sites (variants) we are interested in 
def get_list_of_sites(input_file):
    f = open(input_file)
    head_count = 0  #  For header
    dicti = {}  # Initialize dictionary to keep track of sites
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        # Remove gene name from site id identifier
        full_site_id = data[0]
        site_id_info = full_site_id.split('_')
        site_id = site_id_info[0] + '_' + site_id_info[1] + '_' + site_id_info[2] + '_' + site_id_info[3] + '_' + site_id_info[4]
        dicti[site_id] = 1.0
    return dicti


#  Create dictionary that has keys that are sites, and values that are (nested) dictionaries.
#  The nested dictionaries have keys that are cell_line_ids if that cell_line is heterozygous at the site
#  therefore dictionary can be accessed by is_site_heterozygous[siteID][cellLine].
def learn_heterozygous_sites(sites, het_thresh, het_prob_genotype_file):
    dicti = {}  # Initialize dictionary
    f = open(het_prob_genotype_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#CHROM'):  # If sample id line
            ordered_cell_lines = np.asarray(data[9:])
        if line.startswith('#'):  # Skip other headers
            continue
        site_id = data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3] + '_' + data[4]  # create site_id from line
        if site_id not in sites:  # Ignore sites that we do not have in our data
            continue
        het_probs = np.asarray(data[9:]).astype(float)  # Convert heterozygous probs into array
        # Error checking
        if site_id in dicti:
            print('EROROROR')
            pdb.set_trace()
        dicti[site_id] = {}
        #  Loop through each cellLine
        for i, het_prob in enumerate(het_probs):
            i_cell_line = ordered_cell_lines[i]  # Name of ith cell line
            if het_prob >= het_thresh:  # Is ith cell line heterozygous at this site
                dicti[site_id][i_cell_line] = 1.0
    return dicti


# Main driver for part 3: For various heterozygous probability thresholds, place NA for (sample,site) pairs that have het. prob less than specified threshold
def filter_heterozygous_snps(allelic_counts_gene_mapped_output_file, allelic_counts_gene_mapped_het_only_output_file_root, het_thresh, het_prob_genotype_file):
    #  Create list of sites (variants) we are interested in 
    sites = get_list_of_sites(allelic_counts_gene_mapped_output_file)
    #  Create dictionary that has keys that are sites, and values that are (nested) dictionaries.
    #  The nested dictionaries have keys that are cell_line_ids if that cell_line is heterozygous at the site
    #  therefore dictionary can be accessed by is_site_heterozygous[siteID][cellLine].
    is_site_heterozygous = learn_heterozygous_sites(sites, het_thresh, het_prob_genotype_file)

    # Error checking
    if len(sites) != len(is_site_heterozygous):
        print('ERROR!!! Must stop')
        pdb.set_trace()

    t = open(allelic_counts_gene_mapped_het_only_output_file_root + '.txt', 'w')  # Open output handle
    t_ref = open(allelic_counts_gene_mapped_het_only_output_file_root + '_ref_counts.txt', 'w')  # Open output handle
    t_total = open(allelic_counts_gene_mapped_het_only_output_file_root + '_total_counts.txt', 'w')  # Open output handle
    f = open(allelic_counts_gene_mapped_output_file)  # Open input handle
    head_count = 0  # for header
    # Loop through input file
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # header
            head_count = head_count + 1
            t.write(line + '\n')
            t_ref.write(line + '\n')
            t_total.write(line + '\n')
            #  Make list of ordering of cellLines in allelic_counts files
            ordered_cell_lines = []
            for sample_id in data[1:]:
                cell_line = sample_id.split('_')[0]
                ordered_cell_lines.append(cell_line)
            ordered_cell_lines = np.asarray(ordered_cell_lines)
            continue
        # For each normal line in file
        # Remove gene name from site id identifier
        full_site_id = data[0]
        site_id_info = full_site_id.split('_')
        site_id = site_id_info[0] + '_' + site_id_info[1] + '_' + site_id_info[2] + '_' + site_id_info[3] + '_' + site_id_info[4]

        # ERROR CHECKING
        if site_id not in is_site_heterozygous:
            print('ERRROR: MUST STOP')
            pdb.set_trace()

        # PRINT ONE LINE TO OUTPUT FILE
        t.write(data[0])  # Print row identifier
        # Loop through each sample
        for i,ele in enumerate(data[1:]):
            i_cell_line = ordered_cell_lines[i]  # cell_line of ith sample
            if i_cell_line in is_site_heterozygous[site_id]:  # check if (sample_id, site) is a heterozygous variant at this threshold
                if ele == 'Nan':
                    t.write('\t' + '0_0')
                    t_ref.write('\t' + '0')
                    t_total.write('\t' + '0')
                else:
                    t.write('\t' + ele)
                    t_ref.write('\t' + ele.split('_')[0])
                    t_total.write('\t' + ele.split('_')[1])
            else:
                t.write('\tNan')
                t_ref.write('\tNan')
                t_total.write('\tNan')
        t.write('\n')
        t_ref.write('\n')
        t_total.write('\n')
    f.close()
    t.close()
    t_ref.close()
    t_total.close()


############################################################################################################################################################################################################################
# Input arguments
############################################################################################################################################################################################################################
# Now that we have run the WASP mapping pipeline and GATK ASEReadCounter, we now have one allelic count file per sample
# process_and_organize_allelic_counts.sh will:
######### 1. Merge all of the sample specific count files into one table
######### 2. Map the sites to protein coding genes and remove sites that don't lie on a protein-coding gene
######### 3. For various heterozygous probability thresholds, place NA for (sample,site) pairs that have het. prob less than specified threshold
######### 4. Apply various filters for sites based on number of samples that we have mapped read counts to a het. site (etc)
raw_allelic_counts_dir = sys.argv[1]  # Input dir containing 1 file per sample
processed_allelic_counts_dir = sys.argv[2]  # Ouput dir to save processed matrices
genotype_dir = sys.argv[3]  # Directory containing genotype information
preprocess_total_expression_dir = sys.argv[4]  # Directory containing processed total expression. We use this here because we want to include the sample samples, in the same order as the total expression analysis
gencode_gene_annotation_file = sys.argv[5]  # File containing gencode hg19 gene annotation file. We will use this to map heterozygous sites to genes

############################################################################################################################################################################################################################
############################################################################################################################################################################################################################


# Not all samples are to be used. We performed filtering on samples in the total expression preprocessing.
# Retrieves ordered list of samples to be used (based on total expression filtering)
total_expression_file = preprocess_total_expression_dir + 'quantile_normalized.txt'
ordered_samples = get_ordered_samples_from_total_expression(total_expression_file)

####################################################################################
######### 1. Merge all of the sample specific count files into one table (~10 min)
####################################################################################
raw_allelic_counts_output_file = processed_allelic_counts_dir + 'raw_allelic_counts.txt'  # Output file to save results from step 1
merge_sample_specific_count_files(raw_allelic_counts_dir, ordered_samples, raw_allelic_counts_output_file)


####################################################################################
######### 2. Map the sites to protein coding genes and remove sites that don't lie on a protein-coding gene
####################################################################################
allelic_counts_gene_mapped_output_file = processed_allelic_counts_dir + 'allelic_counts_gene_mapped.txt'  # Output file to save results from step2
map_sites_to_genes(raw_allelic_counts_output_file, allelic_counts_gene_mapped_output_file, gencode_gene_annotation_file)



####################################################################################
######### 3. For various heterozygous probability thresholds, place NA for (sample,site) pairs that have het. prob less than specified threshold
####################################################################################
# Perform for a range of heterozygous thresholds
for het_thresh in [.55, .6, .65, .7, .75, .8, .85, .9, .95,.99,.999]:
    print(het_thresh)
    allelic_counts_gene_mapped_het_only_output_file_root = processed_allelic_counts_dir + 'allelic_counts_gene_mapped_het_prob_' + str(het_thresh).split('.')[1]  # Output root to save results from step3
    het_prob_genotype_file = genotype_dir + 'YRI_het_prob_genotype.vcf'
    filter_heterozygous_snps(allelic_counts_gene_mapped_output_file, allelic_counts_gene_mapped_het_only_output_file_root, het_thresh, het_prob_genotype_file)

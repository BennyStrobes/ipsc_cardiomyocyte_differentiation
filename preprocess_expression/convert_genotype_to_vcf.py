import numpy as np
import os
import sys
import pdb
import gzip


# Create dictionary of all samples we wish to analyze
def get_names_of_samples(sample_names_file): 
    dicti = {}
    f = open(sample_names_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        sample_id = data[0]
        # Cell line for this sample
        cell_line = sample_id.split('_')[0]
        dicti[cell_line] = 1.0
    return dicti

# Main driver function
def convert(input_genotype, vcf_formatted_file, sample_names_file):
    # Create dictionary of all samples we wish to analyze
    sample_dict = get_names_of_samples(sample_names_file)

    # Open input file handle
    f = gzip.open(input_genotype)
    # Open outptu file handle
    t = open(vcf_formatted_file, 'w')

    # Add fileformat tag to top of file
    t.write('##fileformat=VCFv4.2\n')
    # Add date
    t.write('##fileDate=20170817\n')
    # Add format info
    t.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Expected genotype based on imputation">\n')

    # Loop through input file
    for line in f:
        line = line.rstrip()
        data = line.split()

        # Write header without making changes
        if data[0] == '#CHROM':
            # First columns of line should all be kepts
            header = data[0:9]
            # THese columns are sample names
            sample_names = np.asarray(data[9:])
            # Create array of indices of samples that we wish to include
            valid_indices = []
            for i,val in enumerate(sample_names):
                if val in sample_dict:
                    valid_indices.append(i)
            # print line while subsetting to lonly valid sample names
            t.write('\t'.join(header) + '\t' + '\t'.join(sample_names[valid_indices]) + '\n')
            continue

        # The current vcf file has 'nan's for all entries (rows) that are not snvs. We will ignore those
        if data[9] == 'nan':
            continue
        # All other lines..

        # Remove 'chr' prefix to 'chr[num]'
        chrom_num = data[0].split('hr')[1]

        # dosage genotype for each sample
        genos = np.asarray(data[9:])
        # dosage genotypes for the samples we have
        genos_filtered = genos[valid_indices]


        #  Replace value in 'INFO' field (column 8) with '.' (symbolizes missing)
        t.write(chrom_num + '\t' + '\t'.join(data[1:7]) + '\t.\t' + data[8] + '\t' + '\t'.join(genos_filtered) + '\n')

        ####################################################################################
        ####################################################################################
        # VARIOUS CHECKS PERFORMED
        # Make sure 'INFO' wasn't more than just a space filler. It wasn't
        if data[7] != 'INFO':
            print('INFO ASSUMPTION ERROR!')
        # See if every time there was a 'nan', it was for a non-snv. Alwasys held
        if (len(data[3]) != 1 or len(data[4]) != 1) and data[9] != 'nan':
            print('SNV-only ASSUMPTION ERROR!!')

    # Close file handles
    t.close()
    f.close()

#  Helper scipt that extracts dictionary of all sites for genotype analysis
def extract_sites_used(vcf_file):
    dicti = {}
    f = open(vcf_file)
    #  Loop through lines of file
    for line in f:
        line = line.rstrip()
        data = line.split()
        # If line starts with '#', then it is a header and we skip for this purpose
        if line.startswith('#'):
            continue
        # site id is string of form chromNum_position_rsID_refAllele_altAllele
        site_id = data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3] + '_' + data[4]
        #  Add site to dictionary
        dicti[site_id] = 1.0
    return dicti

# Helper method that extracts ordered array of all sample names we are using in our genotype dosage file
def extract_ordered_samples(vcf_file):
    arr = []
    f = open(vcf_file)
    #  Loop through lines of the file
    for line in f:
        line = line.rstrip()
        data = line.split()
        #  This is the one line we care about (so skip if not it)
        if line.startswith('#CHROM') == False:
            continue
        ordered_samples = np.asarray(data[9:])
        return ordered_samples

def convert_het_prob(vcf_formatted_file, vcf_het_prob_file, heterozygous_site_input_dir):
    # Extract dictionary of all sites used in our dosage genotype file
    sites = extract_sites_used(vcf_formatted_file)
    # Extract ordered array of all sample names we are using in our genotype dosage file
    ordered_samples = extract_ordered_samples(vcf_formatted_file)

    # Create array where each element is the index of a sample that we are using for this analysis (used to subset samples)
    f = open(heterozygous_site_input_dir + 'YRI_samples.txt')
    counter = 0
    # Initialize array
    sample_indices = []
    # Loop through samples file
    for i, line in enumerate(f):
        line = line.rstrip()
        data = line.split()
        sample_id = data[0].split('NA')[1]
        if counter >= len(ordered_samples):
            continue
        if ordered_samples[counter] == sample_id:
            sample_indices.append(i)

            counter = counter + 1
    if len(ordered_samples) != len(sample_indices):
        print('SAMPLE ORDERING ERROR (FATAL)')
        pdb.set_trace()
    sample_indices = 3*np.asarray(sample_indices) + 1  # Triplets..
    f.close()


    # Open outptu file handle
    t = open(vcf_het_prob_file, 'w')

    # Add fileformat tag to top of file
    t.write('##fileformat=VCFv4.2\n')
    # Add date
    t.write('##fileDate=20170817\n')
    # Add format info
    t.write('##FORMAT=<ID=HP,Number=1,Type=Float,Description="Heterozygous Probability>\n')

    t.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(ordered_samples) + '\n')

    # Keep track of which sites were used
    used_sites = {}
    # There is a seperate het. prob input file for each chromosome so loop through all chromosomes
    for chrom_num in range(1,23):
        # Het prob file name
        file_name = heterozygous_site_input_dir + 'chr' + str(chrom_num) + '.hg19.impute2.gz'
        f = gzip.open(file_name)
        # Loop through lines (sites) of the file
        for line in f:
            line = line.rstrip()
            data = line.split()
            # create identifier of the site
            # Identifier in same format as  those in sites (chromNum_position_rsID_refAllele_altAllele)
            site_id = str(chrom_num) + '_' + data[2] + '_' + data[1] + '_' + data[3] + '_' + data[4]

            # Remove sites that were filtered out
            if site_id not in sites:
                continue
            used_sites[site_id] = 1.0
            # Write id info to line
            t.write(str(chrom_num) + '\t' + data[2] + '\t' + data[1] + '\t' + data[3] + '\t' + data[4] + '\t100\tPASS\t.\tHP\t')
            # Write het. probs to file
            het_probs = np.asarray(data[5:])[sample_indices]
            t.write('\t'.join(het_probs) + '\n')
        f.close()
    t.close()
    return used_sites

#  Remove sites in vcf_formatted_file that are not in vcf_het_prob_file
def filter_sites(input_file, output_file, sites):
    f = open(input_file)
    t = open(output_file, 'w')
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#'):
            t.write(line + '\n')
            continue
        site_id = data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3] + '_' + data[4]
        if site_id not in sites:
            continue
        t.write(line + '\n')

def convert_het_prob_all_available_samples(vcf_formatted_file, vcf_het_prob_file, heterozygous_site_input_dir):
    # Extract dictionary of all sites used in our dosage genotype file
    sites = extract_sites_used(vcf_formatted_file)
    # Create array where each element is the index of a sample that we are using for this analysis (used to subset samples)
    f = open(heterozygous_site_input_dir + 'YRI_samples.txt')
    counter = 0
    # Initialize array
    sample_indices = []
    ordered_samples = []
    # Loop through samples file
    for i, line in enumerate(f):
        line = line.rstrip()
        data = line.split()
        sample_id = data[0].split('NA')[1]
        ordered_samples.append(sample_id)
        sample_indices.append(i)

    ordered_samples = np.asarray(ordered_samples)
    if len(ordered_samples) != len(sample_indices):
        print('SAMPLE ORDERING ERROR (FATAL)')
        pdb.set_trace()
    sample_indices = 3*np.asarray(sample_indices) + 1  # Triplets..
    f.close()


    # Open outptu file handle
    t = open(vcf_het_prob_file, 'w')

    # Add fileformat tag to top of file
    t.write('##fileformat=VCFv4.2\n')
    # Add date
    t.write('##fileDate=20170817\n')
    # Add format info
    t.write('##FORMAT=<ID=HP,Number=1,Type=Float,Description="Heterozygous Probability>\n')

    t.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(ordered_samples) + '\n')

    # Keep track of which sites were used
    used_sites = {}
    # There is a seperate het. prob input file for each chromosome so loop through all chromosomes
    for chrom_num in range(1,23):
        # Het prob file name
        file_name = heterozygous_site_input_dir + 'chr' + str(chrom_num) + '.hg19.impute2.gz'
        f = gzip.open(file_name)
        # Loop through lines (sites) of the file
        for line in f:
            line = line.rstrip()
            data = line.split()
            # create identifier of the site
            # Identifier in same format as  those in sites (chromNum_position_rsID_refAllele_altAllele)
            site_id = str(chrom_num) + '_' + data[2] + '_' + data[1] + '_' + data[3] + '_' + data[4]

            # Remove sites that were filtered out
            if site_id not in sites:
                continue
            used_sites[site_id] = 1.0
            # Write id info to line
            t.write(str(chrom_num) + '\t' + data[2] + '\t' + data[1] + '\t' + data[3] + '\t' + data[4] + '\t100\tPASS\t.\tHP\t')
            # Write het. probs to file
            het_probs = np.asarray(data[5:])[sample_indices]
            t.write('\t'.join(het_probs) + '\n')
        f.close()
    t.close()
    return used_sites



# Input file provided by Bryce. Close to VCF
input_genotype = sys.argv[1]
# Output_dir
genotype_dir = sys.argv[2]
# Directory containing data relating to herozygous probabilities
heterozygous_site_input_dir = sys.argv[3]
# Names of samples we wish to use
sample_names_file = sys.argv[4]






############################################################################################
#  Convert dosage genotype file to vcf file
############################################################################################


# Output file that contains same information as input_genotype, except in acceptable vcf format
vcf_formatted_file = genotype_dir + 'YRI_genotype_temp.vcf'


# Main driver function
convert(input_genotype, vcf_formatted_file, sample_names_file)




############################################################################################
#  Convert heterozygous probabilities directory to vcf file
############################################################################################

vcf_het_prob_file = genotype_dir + 'YRI_het_prob_genotype.vcf'

# Main driver function to convert heterozygous probs
used_sites = convert_het_prob(vcf_formatted_file, vcf_het_prob_file, heterozygous_site_input_dir)

# Used for debugging sample mis-labelling problem
vcf_het_prob_all_samples_file = genotype_dir + 'YRI_het_prob_genotype_all_samples.vcf'

convert_het_prob_all_available_samples(vcf_formatted_file, vcf_het_prob_all_samples_file, heterozygous_site_input_dir)

############################################################################################
#  Remove sites in vcf_formatted_file that are not in vcf_het_prob_file
#  There are only two sites where this happens.
############################################################################################

vcf_formatted_file_final = genotype_dir + 'YRI_genotype.vcf'
filter_sites(vcf_formatted_file, vcf_formatted_file_final, used_sites)




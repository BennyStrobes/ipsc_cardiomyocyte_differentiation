import numpy as np
import os
import sys
import pdb
import gzip

# Simple helper script to parse genotype_input for only snps on chrom_num_string
# Output is saved to chrom_specific_output_file
def create_chromosome_specific_text_based_snp_file(chrom_num_string, chrom_specific_output_file, genotype_input):
    # Initialize output file handle
    t = open(chrom_specific_output_file, 'w')
    # Initialize input file handle
    f = gzip.open(genotype_input)

    # Loop through lines of input file
    head_count = 0  # used to skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # This line is the header. We can ignore
            head_count = head_count + 1
            continue

        line_chrom_num_string = data[0]  # First element of line is chromosome
        if line_chrom_num_string != chrom_num_string:  # This line contains a SNP on non-desired chromosome --> skip
            continue

        pos = data[1]  # Extract position of the snp
        ref_allele = data[3]  # Extract reference allele of the snp
        alt_allele = data[4]  # Extract alternate allele of the snp

        #  Print pos, ref_allele, alt_allele to space-delimited output file
        t.write(pos + ' ' + ref_allele + ' ' + alt_allele + '\n')

        ##############################################
        # Just some error checking on vcf file from bryce.
        #  Not necessary to run in futre
        ##############################################
        #if data[6] != 'PASS':
        #    print("EROROOROR")
        #    pdb.set_trace()
        #geno = np.asarray(data[9:]).astype(float)
        #maf = sum(geno)/(2.0*len(geno))
        #if maf > .5:
        #    maf = 1 - maf
        #if maf < .01:
        #    print('low MAF')
        #    pdb.set_trace()

    # Close output file handle
    t.close()


# The goal of this script is to create one file for each chromosome that contains only three columns (position, ref_allele, alt_allele)
# These files will be used in WASP mapping pipeline

genotype_input = sys.argv[1]  # Input genotype file
genotype_dir = sys.argv[2]  # Output directory

# Loop through autosomal chromosomes
for chrom_num in range(1, 23):
    print(chrom_num)
    chrom_num_string = 'chr' + str(chrom_num)

    # output file for this specific chromosome
    chrom_specific_output_file = genotype_dir + str(chrom_num) + '.snps.txt'

    # Create this output file.
    create_chromosome_specific_text_based_snp_file(chrom_num_string, chrom_specific_output_file, genotype_input)

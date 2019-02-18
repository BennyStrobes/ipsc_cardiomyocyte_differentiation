import numpy as np 
import os
import sys
import pdb

# Function to extract environmental_variable when environmental_variable_form=='time_steps'
def extract_environmental_variable_time_step_form(sample_name):
    return sample_name.split('_')[1]


# Get Cell line PCs for each sample
def get_cell_line_specific_pc(sample_name, file_name, pc_num):
    cell_liner = sample_name.split('_')[0]
    f = open(file_name)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        if cell_liner == data[0]:
            return data[pc_num]

# Extract sample names from first row of total expression file
def extract_sample_names(total_expression_file):
    f = open(total_expression_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            samples = data[1:]
            continue
    f.close()
    return samples

# Create vectors of names of samples, and their corresponding time step
def get_sample_info(joint_test_input_file):
    # initialize output vectors
    sample_names = []
    time_steps = []
    pc1 = []
    pc2 = []
    pc3 = []
    pc4 = []
    pc5 = []
    pc6 = []
    pc7 = []
    pc8 = []
    pc9 = []
    pc10 = []
    # Stream joint test input file
    f = open(joint_test_input_file)
    head_count = 0  # skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        sample_name = data[0]
        time_step = data[1]
        sample_names.append(sample_name)
        time_steps.append(time_step)
        pc1.append(data[2])
        pc2.append(data[3])
        pc3.append(data[4])
        pc4.append(data[5])
        pc5.append(data[6])
        pc6.append(data[7])
        pc7.append(data[8])
        pc8.append(data[9])
        pc9.append(data[10])
        pc10.append(data[11])
    return np.asarray(sample_names), np.asarray(time_steps), np.asarray(pc1), np.asarray(pc2), np.asarray(pc3), np.asarray(pc4), np.asarray(pc5), np.asarray(pc6), np.asarray(pc7), np.asarray(pc8), np.asarray(pc9), np.asarray(pc10)

# Extract dictionary of test_names (variant_gene pairs), variants, and gene_names
def extract_test_names(target_region_input_file):
    # Initialize output dictionaries
    genes = {}
    variants = {}
    tests = {}
    f = open(target_region_input_file)
    # Stream file
    for line in f:
        line = line.rstrip()
        data = line.split()
        # rs and gene id
        test_info = data[6].split('_')
        rs_id = test_info[1]
        ensamble_id = test_info[2]

        # Add test names to dictionaries
        genes[ensamble_id] = 1
        variants[rs_id] = 1
        tests[rs_id + '_' + ensamble_id] = 1
    return tests, variants, genes


# Now extract gene expression vector for all genes
def extract_gene_expression(genes, total_expression_file, sample_names):
    head_count = 0  # for header
    f = open(total_expression_file)
    for line in f:
        line = line.rstrip()
        data = np.asarray(line.split())
        if head_count == 0:  # header
            # Using header, find ordered indices that correspond to order of sample_names
            head_count = head_count + 1
            reordered_indices = []
            for sample_name in sample_names:
                for i, ele in enumerate(data):
                    if ele == sample_name:
                        reordered_indices.append(i)
            if np.array_equal(data[reordered_indices],sample_names) == False:
                print('ASSUMPTION ERRORRO!')
                pdb.set_trace()
            continue
        # Name of gene for current line
        ensamble_id = data[0]
        # Extract gene expression measurements
        ordered_data = data[reordered_indices].astype(float)
        # Standardize gene expression measurements
        standardized_ordered_data = (ordered_data - np.mean(ordered_data))/(np.std(ordered_data))
        # Check if gene is in our dictionary of genes (used genes)
        if ensamble_id in genes:
            # If it is, add gene expression vector
            genes[ensamble_id] = standardized_ordered_data
    f.close()
    check_to_make_sure_no_missing_entries(genes)
    return genes


def extract_genotype_data(variants, genotype_file, sample_names, genotype_version):
    cell_line_names = sample_to_cell_line_names(sample_names)

    f = open(genotype_file)
    for line in f:
        line = line.rstrip()
        data = np.asarray(line.split())
        if line.startswith('#CHROM'):  # header
            # Using header, find ordered indices that correspond to order of sample_names
            reordered_indices = []
            for cell_line in cell_line_names:
                for i, ele in enumerate(data):
                    if ele == cell_line:
                        reordered_indices.append(i)
            if np.array_equal(data[reordered_indices],cell_line_names) == False:
                print('ASSUMPTION ERRORRO!')
                pdb.set_trace()
            continue
        if line.startswith('#'):
            continue
        # Name of variant for current line
        rs_id = data[2]
        # Get genotype data in same order as sample_names
        ordered_data = data[reordered_indices].astype(float)
        # If we want to round genotypes
        if genotype_version == 'round':
            ordered_data = np.round(ordered_data)
        
        # Check if variant is in our dictionary of variants (used variants)
        if rs_id in variants:
            # If it is, add gene expression vector
            variants[rs_id] = ordered_data
    f.close()
    check_to_make_sure_no_missing_entries(variants)
    return variants

def check_to_make_sure_no_missing_entries(dicti):
    for key in dicti.keys():
        if len(dicti[key]) == 1:
            print('MISSING KEY ERRORR!!')
            pdb.set_trace()
    return

# Convert vector of sample names (cellLine_timeStep) to a vector of cell lines
def sample_to_cell_line_names(sample_names):
    cell_lines = []
    for ele in sample_names:
        cell_lines.append(ele.split('_')[0])
    return np.asarray(cell_lines)

def print_dynamic_eqtl_input_file(output_file, test_names, variants, genes, sample_names, time_steps, environmental_variable_form, pc1, pc2, pc3, pc4, pc5, pc6, pc7,pc8,pc9,pc10):
    # Convert vector of sample names (cellLine_timeStep) to a vector of cell lines
    cell_line_names = sample_to_cell_line_names(sample_names)
    # Open file handle to output file
    t = open(output_file, 'w')
    for test_name in sorted(test_names.keys()):
        # Extract rs_id and ensamble_id from test_name
        rs_id = test_name.split('_')[0]
        ensamble_id = test_name.split('_')[1]
        # extract genotype vector
        geno_vec = variants[rs_id].astype(str)
        # extract total expression vector
        te_vec = genes[ensamble_id].astype(str)
        # Print
        t.write(rs_id + '\t' + ensamble_id + '\t' + ';'.join(time_steps) + '\t' + ';'.join(geno_vec) + '\t' + ';'.join(te_vec) + '\t' + ';'.join(cell_line_names) + '\t' + ';'.join(pc1) + '\t' + ';'.join(pc2) + '\t' + ';'.join(pc3) + '\t' + ';'.join(pc4) + '\t' + ';'.join(pc5) + '\t' + ';'.join(pc6) + '\t' + ';'.join(pc7) + '\t' + ';'.join(pc8) + '\t' + ';'.join(pc9) + '\t' + ';'.join(pc10) + '\n')
    t.close()


##########################
# Command line args
##########################

joint_test_input_file = sys.argv[1]  # Ouput file
dynamic_eqtl_input_file = sys.argv[2]  # Output file
total_expression_file = sys.argv[3]  # File containing processed total expression data
genotype_file = sys.argv[4]  # File with gentoype info for each cell line
environmental_variable_form = sys.argv[5]  # How to encode environmental variable
genotype_version = sys.argv[6]  # Whether to round or use dosage based genotype data
cell_line_specific_pc_file = sys.argv[7]  # File containing PC loadings for each cell line
target_region_input_file = sys.argv[8]  # File containing list of variant gene pairs used in per-time step eqtl analysis. This same set of variant-gene pairs will be used for dynamic eQTL analysis

# Extract sample names from first row of total expression file
sample_names = extract_sample_names(total_expression_file)

# Open output file handle
t = open(joint_test_input_file, 'w')
# Write header to output file handle
t.write('sample_id\tenvironmental_variable\tcell_line_pc1\tcell_line_pc2\tcell_line_pc3\tcell_line_pc4\tcell_line_pc5\tcell_line_pc6\tcell_line_pc7\tcell_line_pc8\tcell_line_pc9\tcell_line_pc10\n')


# Loop through sample names
for sample_name in sample_names:
    # Extract environmental variable (depends on environmental_variable_form)
    if environmental_variable_form == 'time_steps':
        environmental_variable = extract_environmental_variable_time_step_form(sample_name)

    # Get Cell line PCs for each sample
    pc1 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 1)
    pc2 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 2)
    pc3 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 3)
    pc4 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 4)
    pc5 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 5)
    pc6 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 6)
    pc7 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 7)
    pc8 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 8)
    pc9 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 9)
    pc10 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 10)

    # Print information to output file
    t.write(sample_name + '\t' + environmental_variable + '\t' + pc1 + '\t' + pc2 + '\t' + pc3 + '\t' + pc4 + '\t' + pc5 + '\t' + pc6 + '\t' + pc7 + '\t' + pc8 + '\t' + pc9 + '\t' + pc10 + '\n')
t.close()


# Create vectors of names of samples, and their corresponding time step
sample_names, time_steps, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10 = get_sample_info(joint_test_input_file)

# Extract dictionary of test_names (variant_gene pairs), variants, and gene_names
test_names, variants, genes = extract_test_names(target_region_input_file)


# Now extract gene expression vector for all genes
genes = extract_gene_expression(genes, total_expression_file, sample_names)

# Now extract genotype vector for all variants
variants = extract_genotype_data(variants, genotype_file, sample_names, genotype_version)

# Make output file where each row is a test, and columns contain all info (expression, genotype, time steps) to run the test
print_dynamic_eqtl_input_file(dynamic_eqtl_input_file, test_names, variants, genes, sample_names, time_steps, environmental_variable_form, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10)



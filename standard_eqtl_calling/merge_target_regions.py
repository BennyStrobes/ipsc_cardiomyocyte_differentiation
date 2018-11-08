import numpy as np
import os
import sys
import pdb



target_regions_dir = sys.argv[1]
cis_distance = sys.argv[2]
maf_cutoff = float(sys.argv[3])
min_read_count = sys.argv[4]
min_as_read_count = sys.argv[5]
min_het_count = sys.argv[6]


input_root = target_regions_dir + 'target_regions_cis_distance_' + str(cis_distance) + '_maf_cutoff_' + str(maf_cutoff) + '_min_reads_' + str(min_read_count) + '_min_as_reads_' + str(min_as_read_count) + '_min_het_counts_' + str(min_het_count) #+ '_time_' + str(time_step) +'.txt','w')
output_file = input_root + '_merged.txt'

dicti = {}
num_steps = 0
# First loop through all of time steps
# And keep track of the number of time steps a given test passed all of the filters
for time_step in range(0,16):
    num_steps = num_steps + 1
    # Open up target region file for this time step
    input_file = input_root + '_time_' + str(time_step) + '.txt'
    f = open(input_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        namer = data[6]
        # Keep track of in how many time steps each test_name passed the filters
        if namer not in dicti:
            dicti[namer] = 1
        else:
            dicti[namer] = dicti[namer] + 1

# Make new dictionary (all_time_steps) that only has test names found in all time steps
all_time_steps = {}
for key in dicti.keys():
    if dicti[key] == num_steps:
        all_time_steps[key] = 0

# Print the keys in all_time_steps to our output file
t = open(output_file, 'w')
hits = 0
for time_step in range(0,16):
    input_file = input_root + '_time_' + str(time_step) + '.txt'
    f = open(input_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        namer = data[6]
        if namer in all_time_steps and all_time_steps[namer] == 0:
            hits = hits + 1
            all_time_steps[namer] = 1
            t.write(line + '\n')
if hits != len(all_time_steps):
    print('Eroeoror')
import numpy as np
import os
import sys
import pdb




time_step = sys.argv[1]
corrected_quantile_normalized_expression = sys.argv[2]
output_file = sys.argv[3]

f = open(corrected_quantile_normalized_expression)
head_count = 0
for line in f:
    line = line.rstrip()
    data = line.split()
    if head_count == 0:
        head_count = head_count + 1
        sample_names = data[1:]
        continue
t = open(output_file, 'w')
for sample_name in sample_names:
    cell_line = sample_name.split('_')[0]
    curr_time_step = sample_name.split('_')[1]
    if int(time_step) == int(curr_time_step):
        t.write(cell_line + '\n')
t.close()

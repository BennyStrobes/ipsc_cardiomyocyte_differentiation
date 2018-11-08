import numpy as np 
import os
import pdb
import sys
import gzip



def get_samples(genotype_input):
    f = open(genotype_input)
    for line in f:
        line = line.rstrip()
        if line.startswith('#CHROM') == False:
            continue
        data = line.split()
        samples = data[9:]
    return np.asarray(samples)


genotype_input = sys.argv[1]
output_file = sys.argv[2]


samples = get_samples(genotype_input)
t = open(output_file,'w')
for sample in samples:
    t.write(sample + '\n')
t.close()
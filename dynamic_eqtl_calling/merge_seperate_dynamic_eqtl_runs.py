import numpy as np 
import os
import sys
import pdb
import os.path


count = 0

output_file = sys.argv[1]
qtl_results_dir = sys.argv[2]
parameter_string = sys.argv[3]
num_jobs = int(sys.argv[4])



if os.path.isfile(output_file) == True:
    print("Merged output file already exists. Need to delete it before we re-write it")
else:
    t = open(output_file, 'w')

    for job_number in range(num_jobs):
        temp_input_file = qtl_results_dir + parameter_string + '_' + str(job_number) + '.txt'
        f = open(temp_input_file)
        for line in f:
            line = line.rstrip()
            data = line.split()
            if data[-1] == 'NA':
                count = count + 1
                continue
            t.write(line + '\n')
        f.close()
        os.system('rm ' + temp_input_file)
    print(str(count) + ' tests have NAs')
    t.close()

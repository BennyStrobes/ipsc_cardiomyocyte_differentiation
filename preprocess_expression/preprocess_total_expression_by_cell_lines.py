import numpy as np 
import os
import sys
import pdb





def get_unique_cell_lines(quantile_normalized_data):
    f = open(quantile_normalized_data)
    head_count = 0
    cell_lines = []
    time_steps = []
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            for ele in data[1:]:
                cell_line = ele.split('_')[0]
                time_step = ele.split('_')[1]
                cell_lines.append(cell_line)
                time_steps.append(time_step)
            continue
    return np.unique(cell_lines), np.asarray(cell_lines), np.asarray(time_steps)




# Print cell line expression to output file
def create_cell_line_expression_ignore_missing(unique_cell_lines, ordered_cell_lines, ordered_time_steps, quantile_normalized_data, cell_line_expr_ignore_missing_file):
    # create mapping from cell line name to position
    cell_line_mapping = {}
    for i, cell_line in enumerate(unique_cell_lines):
        cell_line_mapping[cell_line] = i

    # open input file
    f = open(quantile_normalized_data)
    # open output file
    t = open(cell_line_expr_ignore_missing_file, 'w')
    # Print header
    t.write('gene_name\t' + '\t'.join(unique_cell_lines) + '\n')
    head_count = 0
    # Stream input file
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[0]
        expr = np.asarray(data[1:]).astype(float)
        for time_step in range(16):
            time_step_specific_expr = np.zeros(len(unique_cell_lines))
            counter = 0
            time_step_str = str(time_step)
            for pos, expr_val in enumerate(expr):
                curr_time = ordered_time_steps[pos]
                curr_line = ordered_cell_lines[pos]
                if curr_time != time_step_str:
                    continue
                cell_line_pos = cell_line_mapping[curr_line]
                time_step_specific_expr[cell_line_pos] = expr_val
                counter = counter + 1
            if counter != len(unique_cell_lines):
                continue
            standardized = (time_step_specific_expr - np.mean(time_step_specific_expr))/np.std(time_step_specific_expr)
            t.write(ensamble_id + '_' + time_step_str + '\t' + '\t'.join(standardized.astype(str)) + '\n')
    t.close()




preprocess_total_expression_dir = sys.argv[1]

quantile_normalized_data = preprocess_total_expression_dir + 'quantile_normalized_no_projection.txt'

unique_cell_lines, ordered_cell_lines, ordered_time_steps = get_unique_cell_lines(quantile_normalized_data)


# Output file
cell_line_expr_ignore_missing_file = preprocess_total_expression_dir + 'cell_line_expression_ignore_missing.txt'
create_cell_line_expression_ignore_missing(unique_cell_lines, ordered_cell_lines, ordered_time_steps, quantile_normalized_data, cell_line_expr_ignore_missing_file)


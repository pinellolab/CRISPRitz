#!/usr/bin/env python

'''
From a clustered file, extract the first x Top1.
This version of the script is DIFFERENT to the CRISPRme version (add #!/usr/bin/env python to first line, print the max mm, DNA and RNA bulges to be used in summary by guide pos script )
'''

#Script extract only the top X from all clusters
# argv1 is input file ordered in cluster
# argv2 is job_id 
import sys
header = {'Bulge Type':0, '#Bulge Type':0, 'Bulge type':0,
        'crRNA':1, 'DNA':2,
        'Chromosome':3, 'Position':4, 'Cluster Position':5 ,'Direction':6,
        'Mismatches':7, 'Bulge Size':8, 'Total':9,'Min_mismatches':10, 'Max_mismatches':11,
        'PAM_disr':12, 'PAM_gen':13, 'Var_uniq':14, 'Samples':15}

top_x = 1
top_x_inserted = 0
job_id = sys.argv[2]
max_mm = 0
max_dna = 0
max_rna = 0
with open(sys.argv[1]) as targets, open( job_id + '.top_' + str(top_x) + '.txt', 'w+') as result:
    line = targets.readline()
    if '#' in line:
        result.write(line)  #header
        line = targets.readline().strip().split('\t')
    else:
        line = line.strip().split('\t')
    prev_position = line[header['Cluster Position']]
    prev_chr = line[header['Chromosome']]
    prev_direction = line[header['Direction']]
    max_mm = int(line[header['Mismatches']])
    if line[0] == 'DNA':
        max_dna = int(line[header['Bulge Size']])
    elif line[0] == 'RNA':
        max_rna = int(line[header['Bulge Size']])

    result.write('\t'.join(line) + '\n')
    top_x_inserted = top_x_inserted + 1
    for line in targets:
        line = line.strip().split('\t')
        current_pos = line[header['Cluster Position']]
        current_chr = line[header['Chromosome']]
        current_direction = line[header['Direction']]
        if int (line[header['Mismatches']]) > max_mm:
            max_mm = int(line[header['Mismatches']])
        if line[0] == 'DNA':
            if int(line[header['Bulge Size']]) > max_dna:
                max_dna = int(line[header['Bulge Size']])
        elif line[0] == 'RNA':
            if int(line[header['Bulge Size']]) > max_rna:
                max_rna = int(line[header['Bulge Size']])


        #if i can still instert top1
        if prev_position == current_pos and prev_chr == current_chr and prev_direction == current_direction and top_x_inserted < top_x:
            top_x_inserted = top_x_inserted + 1
            result.write('\t'.join(line) + '\n')
        elif prev_position != current_pos or prev_chr != current_chr or prev_direction != current_direction:
            top_x_inserted = 1
            prev_position = current_pos
            prev_chr = current_chr
            prev_direction = current_direction
            result.write('\t'.join(line) + '\n')
# print(str(max_mm) + ',' + str(max_dna) + ',' + str(max_rna))
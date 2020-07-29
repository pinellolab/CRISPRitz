#!/usr/bin/env python


'''
Script for creating dictionaries for later sample analysis. The name of the fasta file and the name in the #CHROM column MUST be the same.
Use http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ as vcf reference
WARNING! About 20 minutes per vcf file, for a total of 80GB of disk space
'''

import gzip
import sys
import json
import time

# argv1 is ALLchrx.gzip --> Create the json file containing chr_pos -> 'list,of,sample;reference_char,var_char'
# argv2 is chr name --> chr21
# argv3 is directory to store the dictinary
chr1_dict = dict()
start_time = time.time()
add_to_name = ''     #string to add to the chr number, eg in VCF hg38 -> is already chr1, so add_to_name is '';
                        #VCF hg19 -> is 1, so add_to_name is 'chr'

with gzip.open(sys.argv[1], 'rb') as targets:
    for line in targets:
        line = line.decode('ascii')
        if ('#CHROM') in line:
            column_vcf = line.strip().split('\t')
            break
    for line in targets:                #Save CHROM [0], POS[1], REF [3], ALT [4], List of Samples [9:]
        line = line.decode('ascii').strip().split('\t')
        list_samples = []
        list_chars = []
        for pos, i in enumerate(line[9:]):
            if ('1' in i):
                list_samples.append(column_vcf[ pos + 9])
        string_test = add_to_name + line[0] + ',' + line[1]
        #Add in last two position the ref and alt nucleotide, eg: chrX,100 -> sample1,sample5,sample10;A,T
        #If no sample was found, the dict is chrX,100 -> ;A,T
        list_chars.append(line[3])
        list_chars.append(line[4])
        try:
            chr1_dict[string_test] = ','.join(list_samples) + ';' + ','.join(list_chars)
        except:
            chr1_dict[string_test] = ';' + ','.join(list_chars) #None
        #result.write(line[0] + '\t' + line[1] + '\t' + line[3] + '\t' + line[4] + '\t' + ','.join(list_samples) + '\n')
        
with open(sys.argv[3] + '/my_dict_' + sys.argv[2] + '.json', 'w') as f:
    json.dump(chr1_dict, f) 
print('Created ' + 'my_dict_' + sys.argv[2] + '.json' + ' in', time.time() - start_time)
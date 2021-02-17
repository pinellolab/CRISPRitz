#!/usr/bin/env python

'''
Script for summary by guide. Takes in input the samples.annotated.txt and count the targets occurencies for every mm/bulge value.
Crea anche file per tabella generale guide: se un target ha sample lo conto come enriched. 
Se due target sono dello stesso gruppo scomposto (chr pos uguale di seguito), per la tabella generale conto solo il primo, per summary
by guide li conto entrambi
'''

# sys1 samples.annotation.txt
# sys2 mms
# sys3 bulges DNA
# sys4 bulges RNA
# sys5 guide file
# sys6 is jobid
# sys7 is type of post-process done ('No' -> no post process done, cannot count uniq_var | 'Uniq' -> post process done, can count uniq_var)

import sys
import numpy as np
import pandas as pd
import subprocess


mms = int(sys.argv[2])
bulges_dna = int(sys.argv[3])
bulges_rna = int(sys.argv[4])
guide_file = sys.argv[5]
type_post = sys.argv[7]
bulge = bulges_dna
if bulges_rna > bulges_dna:
    bulge = bulges_rna
guide_dict = dict()
count_creation = dict() #{GUIDE1 -> [0 0 0;0 0 0] per ogni categoria mms-bulge,}

general_table = dict() 
'''
{GUIDE1 -> 
    {
        ref -> [0 0 0 0 0 0],      #len = mms + max_bulge     
        var -> [0 0 0 0 0 0],      #len = mms + max_bulge
    }
}
'''
with open (guide_file,'r') as all_guides:
    for line in all_guides:
        line = line.strip()
        # if type_post == 'Uniq':
        guide_dict[line] = [np.zeros((mms + 1, 1)), np.zeros((mms + 1,bulges_dna + 1)), np.zeros((mms + 1,bulges_rna + 1)),
                            np.zeros((mms + 1, 1)), np.zeros((mms + 1,bulges_dna + 1)), np.zeros((mms + 1,bulges_rna + 1))]
        count_creation[line] = [np.zeros((mms + 1, 1)), np.zeros((mms + 1,bulges_dna + 1)), np.zeros((mms + 1,bulges_rna + 1))]
        general_table[line] = dict()
        general_table[line]['ref'] = np.zeros(mms + bulge + 1, dtype = int)
        general_table[line]['var'] = np.zeros(mms + bulge + 1, dtype = int)
        # else:   #TODO aggiungere nel caso di ricerca solo var
        #     guide_dict[line] = [np.zeros((mms + 1, 1)), np.zeros((mms + 1,bulges_dna + 1)), np.zeros((mms + 1,bulges_rna + 1))]
        #     general_table[line] = dict()
        #     general_table[line]['ref'] = np.zeros(mms + bulge + 1, dtype = int)
        #     general_table[line]['var'] = np.zeros(mms + bulge + 1, dtype = int)
current_scomposition = ''
with open(sys.argv[1]) as targets:    
         
    for line in targets:            
        if '#' in line:
            continue
        line = line.strip().split('\t')
        
        if line[0] == 'X':
            if line[12] != 'n':     #If there are samples, add to enriched count, else add to reference count (SUM BY GUIDE)
                guide_dict[line[1].replace('-','')][3][int(line[7])][int(line[8])] += 1
            else:
                guide_dict[line[1].replace('-','')][0][int(line[7])][int(line[8])] += 1 
        elif line[0] == 'DNA':
            if line[12] != 'n':
                guide_dict[line[1].replace('-','')][4][int(line[7])][int(line[8])] += 1
            else:
                guide_dict[line[1].replace('-','')][1][int(line[7])][int(line[8])] += 1
        else:
            if line[12] != 'n':
                guide_dict[line[1].replace('-','')][5][int(line[7])][int(line[8])] += 1
            else:
                guide_dict[line[1].replace('-','')][2][int(line[7])][int(line[8])] += 1
        if line[1].replace('-','') + line[3] + line[5] + line[6] != current_scomposition:
            #New scomposition
            current_scomposition = line[1].replace('-','') + line[3] + line[5] + line[6]   
            if line[12] != 'n':         #sample associated, meaning it's an enriched target (GENERAL TABLE)
                general_table[line[1].replace('-','')]['var'][int(line[7]) + int(line[8])] += 1
            else:
                general_table[line[1].replace('-','')]['ref'][int(line[7]) + int(line[8])] += 1

        #Count pam creation: -> Currently count is done on every line #Count only on first of scomposition
        if line[10] != 'n':
            if line[0] == 'X':
                count_creation[line[1].replace('-','')][0][int(line[7])][int(line[8])] += 1 
            elif line[0] == 'DNA':
                count_creation[line[1].replace('-','')][1][int(line[7])][int(line[8])] += 1 
            else:
                count_creation[line[1].replace('-','')][2][int(line[7])][int(line[8])] += 1 
    
    if type_post == 'Uniq':
        #Open addToGeneralTable.txt and add to each guide the values of 'ref' semicommon
        with open (sys.argv[6] + '.addToGeneralTable.txt') as add_values:
            for ref_semicommon_counts in add_values:
                ref_semicommon_counts = ref_semicommon_counts.strip().split('\t')   #GUIDE 0total 1total 2total ... 10total
                for i in range(mms + bulge + 1):
                    general_table[ref_semicommon_counts[0]]['ref'][i] = general_table[ref_semicommon_counts[0]]['ref'][i] + int(ref_semicommon_counts[i + 1]) 

    #Save results
    with open(sys.argv[6] + '.general_target_count.txt', 'w+') as general_count:
        count_for = '(' + ' - '.join([str(tot) for tot in range (1, mms + bulge + 1)]) + ' Mismatches + Bulges)'
        general_count.write('#Guide\tOn-Targets (Reference - Enriched)\tOff-Targets Reference ' + count_for + '\tOff-Targets Enriched ' + count_for + '\n')
        for guide in guide_dict.keys():
            tab_summary = pd.DataFrame(columns = ['Guide', 'Bulge Type', 'Bulge Size', 'Mismatches', 'Targets in Reference', 'Targets in Enriched', 'PAM Creation'])
            for m in range(mms + 1):
                for b_d in range(bulges_dna +1):
                    tab_summary =tab_summary.append({'Guide': guide, 'Bulge Type': 'DNA', 'Bulge Size': b_d, 'Mismatches': m, 'Targets in Reference': guide_dict[guide][1][m][b_d], 'Targets in Enriched':guide_dict[guide][4][m][b_d], 'PAM Creation':count_creation[guide][1][m][b_d]  }, ignore_index = True)            

                for b_r in range(bulges_rna +1):
                    tab_summary =tab_summary.append({'Guide': guide, 'Bulge Type': 'RNA', 'Bulge Size': b_r, 'Mismatches': m, 'Targets in Reference': guide_dict[guide][2][m][b_r], 'Targets in Enriched': guide_dict[guide][5][m][b_r], 'PAM Creation':count_creation[guide][2][m][b_r] }, ignore_index = True)

                tab_summary =tab_summary.append({'Guide': guide, 'Bulge Type': 'X', 'Bulge Size': 0, 'Mismatches': m, 'Targets in Reference': guide_dict[guide][0][m][0], 'Targets in Enriched':guide_dict[guide][3][m][0], 'PAM Creation':count_creation[guide][0][m][0] }, ignore_index = True)
            tab_summary.to_csv(sys.argv[6] + '.summary_by_guide.' + guide +'.txt', index = False, sep = '\t')
            
            general_count.write(guide + '\t' + str(general_table[guide]['ref'][0] + general_table[guide]['var'][0]) + ' (' + str(general_table[guide]['ref'][0]) + ' - ' + str(general_table[guide]['var'][0]) + ')\t' +
                                str(sum(general_table[guide]['ref'][1:])) + ' (' + ' - '.join([ str(x) for x in general_table[guide]['ref'][1:]]) + ')\t' + 
                                str(sum(general_table[guide]['var'][1:])) + ' (' + ' - '.join([ str(x) for x in general_table[guide]['var'][1:]]) + ')\n'                             
                                )

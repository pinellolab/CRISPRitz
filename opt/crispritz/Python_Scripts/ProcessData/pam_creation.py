#!/usr/bin/env python

'''
Create the PAM Creation column
This version of the script is DIFFERENT to the CRISPRme version (argv4 is output name, modified ext_seq_pam_creation call, modified bedfile name)
'''

# sys argv 1 is input file (uniq.sorted.txt)
# sys argv 2 is pam NOTE tested only with NRG
# sys argv 3 is reference genome dir
# sys argv 4 is output name
# First scan of file: when IUPAC in pam, save position into bedfile
# Use bedtools getfasta to extract all ref sequences
# Second scan of file: add column with 
# -; G; -;
# if from C to ref (and it was not a pam), we got S in enr (created new pam).
# Not sure if calling bedtools for each iupac pam found + single scan of file is faster than single bedtools call + 2 scan of file 

import sys
import subprocess
import os 
from os import listdir                      #for getting directories
from os.path import isfile, isdir,join      #for getting directories
import itertools
import time
pam = sys.argv[2].strip().upper()
guide_len = 20
pam_at_beginning = False
genome_dir = os.path.realpath(sys.argv[3])
with open (sys.argv[2]) as pam:
    line = pam.read().strip()
    pam = line.split(' ')[0]
    len_pam = int(line.split(' ')[1])
    guide_len = len(pam) - len_pam
    pos_beg = 0
    pos_end = None
    pam_begin = 0
    pam_end = len_pam * (-1)
    if len_pam < 0:
        guide_len = len(pam) + len_pam
        pam = pam[: (len_pam * (-1))]
        len_pam = len_pam * (-1)
        pos_beg = len_pam
        pos_end = None
        pam_begin = 0
        pam_end = len_pam
        pam_at_beginning = True
    else:
        pam = pam[(len_pam * (-1)):]
        pos_beg = 0
        pos_end = len_pam * (-1)
        pam_begin = len_pam * (-1)
        pam_end = None

name_output = sys.argv[4]#sys.argv[1][:sys.argv[1].rfind('.')]

iupac_code = {
          "R":("A", "G"),
          "Y":("C", "T"),
          "S":("G", "C"),
          "W":("A", "T"),
          "K":("G", "T"),
          "M":("A", "C"),
          "B":("C", "G", "T"),
          "D":("A", "G", "T"),
          "H":("A", "C", "T"),
          "V":("A", "C", "G"),
          "r":("A", "G"),
          "y":("C", "T"),
          "s":("G", "C"),
          "w":("A", "T"),
          "k":("G", "T"),
          "m":("A", "C"),
          "b":("C", "G", "T"),
          "d":("A", "G", "T"),
          "h":("A", "C", "T"),
          "v":("A", "C", "G"),
          'N':('A', 'T', 'C', 'G')
          }

iupac_code_set = {
          "R":{"A", "G"},
          "Y":{"C", "T"},
          "S":{"G", "C"},
          "W":{"A", "T"},
          "K":{"G", "T"},
          "M":{"A", "C"},
          "B":{"C", "G", "T"},
          "D":{"A", "G", "T"},
          "H":{"A", "C", "T"},
          "V":{"A", "C", "G"},
          "r":{"A", "G"},
          "y":{"C", "T"},
          "s":{"G", "C"},
          "w":{"A", "T"},
          "k":{"G", "T"},
          "m":{"A", "C"},
          "b":{"C", "G", "T"},
          "d":{"A", "G", "T"},
          "h":{"A", "C", "T"},
          "v":{"A", "C", "G"},
          "A":{"A"},
          "T":{"T"},
          "C":{"C"},
          "G":{"G"},
          "a":{"a"},
          "t":{"t"},
          "c":{"c"},
          "g":{"g"},
          'N':{'A','T','G','C'}
        }
rev_comp_pam = dict()
for i in range(len_pam):
    rev_comp_pam[i] = len_pam - i - 1
# rev_comp_pam = {
#     0 : 2,
#     1: 1,
#     2: 0
# }

def rev_comp(a):
    if a == 'A' or a == 'a':
        return 'T'
    if a == 'T' or a == 't':
        return 'A'
    if a == 'C' or a == 'c':
        return 'G'
    return 'C'
start_time = time.time()
#Bedtools estrae da start(compreso) a end(non compreso)
if not pam_at_beginning:
    with open (sys.argv[1]) as uniq, open(name_output + '.bed', 'w+') as bedfile:
        #header=uniq.readline()   #NOTE uncomment if file has header
        for line in uniq:
            if '#' in line:
                continue
            line = line.strip().split('\t')
            for pos, char in enumerate(line[2][pam_begin:pam_end]):           
                if char in iupac_code:
                    if line[0] == 'DNA':
                        bulges = int(line[8])
                    elif line[0] == 'RNA':
                        bulges = int(line[8]) * (-1)       
                    else:
                        bulges = 0
                    if line[6] == '+':
                        # La pos (line[4] mi indica la prima lettera del target, eg la prima A in ACGCGACTAGCTACGCACGTNRG)
                        bedfile.write(line[3] + '\t' + str( int(line[4]) + guide_len + bulges + pos) + '\t' + str(int(line[4]) + guide_len + bulges + pos + 1) +'\n') 
                    else:
                        # La pos (line[4] mi indica la prima lettera del target, che sarebbe la prima lettera della pam rev compl, 
                        # eg la prima C in CCGTGCATACTAGCTACGCACGT)
                        bedfile.write(line[3] + '\t' + str(int(line[4]) + rev_comp_pam[pos]) + '\t' + str(int(line[4]) + rev_comp_pam[pos] + 1) +'\n')
else:
    with open (sys.argv[1]) as uniq, open(name_output + '.bed', 'w+') as bedfile:
        #header=uniq.readline()   #NOTE uncomment if file has header
        for line in uniq:
            if '#' in line:
                continue
            line = line.strip().split('\t')
            for pos, char in enumerate(line[2][pam_begin:pam_end]):           
                if char in iupac_code:
                    if line[0] == 'DNA':
                        bulges = int(line[8])
                    elif line[0] == 'RNA':
                        bulges = int(line[8]) * (-1)        
                    else:
                        bulges = 0
                    if line[6] == '+':
                        #SEE file pamatbeginning.txt
                        bedfile.write(line[3] + '\t' + str( int(line[4]) + pos) + '\t' + str(int(line[4]) + pos + 1) +'\n') 
                    else:
                        # #SEE file pamatbeginning.txt
                        bedfile.write(line[3] + '\t' + str(int(line[4]) + guide_len + bulges + rev_comp_pam[pos]) + '\t' + str(int(line[4]) + guide_len + bulges + rev_comp_pam[pos] + 1) +'\n')

print('Created bedfile:', time.time() - start_time)
start_time = time.time()
bed_file = os.path.realpath(name_output + '.bed')
subprocess.call(['./ext_seq_pam_creation.sh ' + genome_dir + ' ' + bed_file], shell = True) 
print('Bedtools DONE:', time.time() - start_time)
onlyfile = [f for f in listdir('tmp_seq') if isfile(join('tmp_seq', f))]
start_time = time.time()
open_file = dict()

for i in onlyfile:
    open_file[i.split('.fa')[0]] =  open('tmp_seq/' + i)
    #open_file.append([i, open('tmp_seq/' + onlyfile)])

# for i in range (1,23):
#     try:
#         open_file.append(open('chr' + str(i) + '.fa_seq.txt'))
#     except:
#         open_file.append(open('chr' + '1' + '.fa_seq.txt')) #TODO sistemare per indicare che se non ho iupac nella pam non ho questo file
# try:
#     open_file.append(open('chrX.fa_seq.txt'))
# except:
#     open_file.append(open('chr' + '1' + '.fa_seq.txt')) #TODO sistemare per indicare che se non ho iupac nella pam non ho questo file
# count21 = 0
# count22 = 0
with open(sys.argv[1]) as uniq, open(name_output + '.pamcreation.txt', 'w+') as res:
    for line in uniq:
        if '#' in line:
            continue
        line = line.strip().split('\t')
        file_corresponding = line[3]
        # if file_corresponding == 'chr21':
        #     count21 = count21 + 1
        # elif file_corresponding == 'chr22':
        #     count22 = count22 + 1
        total_line = []
        found_creation = False
        # if int(line[4]) == 22687078:
        #     print(line[2])
        for pos, char in enumerate(line[2][pam_begin:pam_end]):  
            # if int(line[4]) == 22687078:
            #     print('Char:', char)
            if char in iupac_code:
                ref_char = open_file[file_corresponding].readline().strip().upper()
                # print('file_corresponding', file_corresponding)
                # print('ref_char', ref_char)
                # print('21', count21)
                # print('22', count22)
                if line[6] == '-':
                    ref_char = rev_comp(ref_char)
                # if int(line[4]) == 22687078:
                #     print('Ref char:' , ref_char)
                
                char_to_write = iupac_code_set[pam[pos]] & iupac_code_set[ref_char]
                
                # if int(line[4]) == 22687078:
                #     print('Char to write:', char_to_write)
                if not char_to_write:
                    #char_to_write =   iupac_code_set[char] - iupac_code_set[ref_char] #char_to_write
                    # if int(line[4]) == 22687078:
                    #     print('Char to write2: ', char_to_write)
                    found_creation = True
                char_to_write = iupac_code_set[pam[pos]] & iupac_code_set[char]
                #char_to_write = ','.join(char_to_write)
                total_line.append(char_to_write)
            else: 
                total_line.append(char)
            # if int(line[4]) == 22687078:
            #     print('Total_line', total_line)
        
        #Total
        #line.append(str(int(line[6]) + int(line[7])))
         #Min
        line.append('-')
        #max
        line.append('-')
        #PAM dis
        line.append('n')
        #Pam create
        if found_creation:
            list_pam = []
            for el in itertools.product(*total_line):
                list_pam.append(''.join(el))
            line.append(','.join(list_pam))
        else:
            line.append('n')
        #Var uniq
        line.append('y')
        
        res.write('\t'.join(line) + '\n') 



# with open(sys.argv[1]) as uniq, open(name_output + '.pamcreation.txt', 'w+') as res:
#     for line in uniq:
#         line = line.strip().split('\t')
#         file_corresponding = line[3].split('chr')[-1]
#         if file_corresponding == 'X':
#             file_corresponding = -1
#         else:
#             file_corresponding = int(file_corresponding) -1
#         total_line = ''
#         for pos, char in enumerate(line[2][-3:]):           #TODO depend on PAM size, begin/end
#             if char in iupac_code:
#                 ref_char = open_file[file_corresponding].readline().strip().upper()
#                 if line[5] == '-':
#                     ref_char = rev_comp(ref_char)
#                 char_to_write =   iupac_code_set[char] - iupac_code_set[ref_char]
#                 char_to_write = ','.join(char_to_write)
#                 total_line = total_line + char_to_write +';'
#             else: 
#                 total_line = total_line + '-;'
            
#         line.append(str(int(line[6]) + int(line[7])))
#          #Min
#         line.append('-')
#         #max
#         line.append('-')
#         #PAM dis
#         line.append('n')
#         #Pam create
#         line.append(total_line)
#         #Var uniq
#         line.append('y')
        
#         res.write('\t'.join(line) + '\n') 

# for i in range(len(open_file)):
#     open_file[i].close

print('Done', time.time() - start_time)

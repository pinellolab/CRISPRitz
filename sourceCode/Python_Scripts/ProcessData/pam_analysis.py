#!/usr/bin/env python



'''
Add Min Max PAM Disruption columns
This version of the script is DIFFERENT to the CRISPRme version (argv4 is output name)
'''

#Was called min_max.py
#Calculates, for the semicommon file, the total column, the min max mismatches column, the pam disr column, pam creation (with 'n' value), var uniq (with 'n' value)
#sys 1 is semicommon file
#sys 2 is pam file NOTE not tested with different from NRG

# If a iupac char is lowercase, it means that every real value it corresponds is a mismatch, meaning that changing it with a real char does
# not add a new theoretical mismatch (eg:
# guide:    AATCCTAG...
# target:   AArCCTAG..
# targetA:  AAaCCTAG...
# targetG:  AAgCCTAG...
# )
# If a iupac char is uppercase, it means that at least one of the corresponding real char is a mismatch (NOTE true only if the guide has no iupac
# char), meaning that there is a theoretical new mismatch (eg:
# guide:    AATCCTAG...
# target:   AAYCCTAG..
# targetC:  AAcCCTAG... > new mismatch
# targetT:  AATCCTAG...
# )
# TODO nel caso della pam, se ho uno iup nel targ e uno nella guida, prendo quello che non rimane nell'intersezione dei loro valori nel dict
# se ho iup nel tar e normale nella guida, prendo i valori diversi dalla guida nei valori del dict (NOTE faccio set(target) - set(guida))
# eg    NRG
#       NMG     -> segnalo C
#       NGG
#       NBG     -> segnalo C,T
import sys
import itertools
#argv 1 is cluster file
#argv 2 is pam file
#argv 3 is type of genome: var or both
#argv4 is output name
fill_column = True
if sys.argv[3] == 'var':
    fill_column = False
with open (sys.argv[2]) as pam:
    line = pam.read().strip()
    pam = line.split(' ')[0]
    len_pam = int(line.split(' ')[1])
    pos_beg = 0
    pos_end = None
    pam_begin = 0
    pam_end = len_pam * (-1)
    if len_pam < 0:
        pam = pam[: (len_pam * (-1))]
        len_pam = len_pam * (-1)
        pos_beg = len_pam
        pos_end = None
        pam_begin = 0
        pam_end = len_pam
    else:
        pam = pam[(len_pam * (-1)):]
        pos_beg = 0
        pos_end = len_pam * (-1)
        pam_begin = len_pam * (-1)
        pam_end = None

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

name_output = sys.argv[4] #sys.argv[1][:sys.argv[1].rfind('.')]
with open (sys.argv[1]) as t, open(name_output + '.minmaxdisr.txt','w+') as result:
    header = '#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge Size\tTotal\tMin_mismatches\tMax_mismatches\tPam_disr'
    if fill_column:
        header += '\tPAM_gen\tVar_uniq'
    header += '\n'
    result.write(header)
    next(t) #Skip header
    for line in t:
        found_iupac = False
        found_iupac_pam = False
        line = line.strip().split('\t')
        max_mm = int(line[7])
        #for pos, char in enumerate(line[2]):
        
        for char in line[2][pos_beg: pos_end]:      
            if char in iupac_code:      
                found_iupac = True
                # if any(line[1][pos] != x for x in iupac_code[char]):        #TODO add IUPAC guide char support: a=line[1][pos] -> for each chr in iupac_code[a]
                #     max_mm = max_mm + 1  
                if char.isupper():
                    max_mm = max_mm + 1
    
        #Pam disruption 
        pam_disr = []
        for pos, char in enumerate(line[2][pam_begin:pam_end]):
            if char in iupac_code:      
                
                diff_res = iupac_code_set[char] - iupac_code_set[pam[pos]]
                if diff_res:
                    pam_disr.append(diff_res)
                    found_iupac_pam = True
                else:
                    pam_disr.append(iupac_code_set[char])
                

                # if pam[pos] in iupac_code:      #both are IUPAC
                #     for normal_pam_char in iupac_code[pam[pos]]:
                #         for normal_target_char in iupac_code[char]:
                #             if normal_pam_char != normal_target_char:
                #                 corr_targ_disr = corr_targ_disr + normal_target_char + ','
                #         if corr_targ_disr:    
                #             pam_disr = pam_disr +  '(' + str(pos) + ',' + normal_pam_char + ') -> (' + corr_targ_disr[:-1] + ')' + ';'
                #         corr_targ_disr = ''
                # else: #confronto normale
                #     for normal_target_char in iupac_code[char]:
                #         if pam[pos] != normal_target_char:
                #             corr_targ_disr = corr_targ_disr + pam[pos] + ','
                #     if corr_targ_disr:
                #        pam_disr = pam_disr +  '(' + str(pos) + ',' + pam[pos] + ') -> (' + corr_targ_disr[:-1] + ')' + ';' 
                #     corr_targ_disr = ''
            else:
                pam_disr.append(char)
        #print (pam_disr)
        if found_iupac:
            line.append(line[7])
            line.append(str(max_mm))
        else:
            line.append('-')
            line.append('-')
        if (found_iupac_pam):
            pam_disr_list = []
            for p in itertools.product(*pam_disr):
                pam_disr_list.append(''.join(p))
            line.append(','.join(pam_disr_list))
        else:
            line.append('n')
        #Pam creation
        if fill_column:
            line.append('n')
        #Uniq var
            line.append('n')
        result.write('\t'.join(line) + '\n')
        
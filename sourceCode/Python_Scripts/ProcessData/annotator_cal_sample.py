#!/usr/bin/env python

'''
Merge of annotator, calc_samples_faster.py and scores
Prende in input il file dei top1, ordinati per chr, e estrae i samples corrispondenti. Per ogni target, salva l'insieme dei sample in samples.all.txt, crea le combinazioni tenendo i target reali
in samples.txt, poi calcola l'annotazione corrispondente e crea il file Annotation.targets e  i vari summaries. . Salva in samples.annotation I REF, i VAR degli Uniq, i VAR e i REF dei Semicommon
Added compatibility with dictionary chr_pos -> s1,s2;A,C/sNew;A,T
'''


#NOTE serve 20130606_sample_info.xlsx nella stessa cartella di questo script 
#argv1 è il file .bed con le annotazioni
#argv2 è il file .cluster.txt, che è ordinato per cromosoma. Era (03/03) il file top1 ordinato per chr
#argv3 è nome del file in output
#argv4 è directory dei dizionari
#argv5 is pamfile
#argv 6 is max allowed mms
#argv 7 is genome reference directory (Eg ../../Genomes/hg38_ref)
#argv8 is guide file
#argv9 is max allowed DNA bulges
#argv10 is max allowed RNA bulges
#argv11 is sampleID file for sample - population - superpopulation association
#argv -1 is step
# NOTE 06/03  -> removed PAM Disruption calculation
#NOTE 29/03 -> le colonne min max sono rimosse, dal file total.cluster sono già presenti colonne sample, annotation, real guide
# 29/03 colonne in input #Bulge_type     crRNA   DNA     Chromosome      Position        Cluster Position        Direction       Mismatches      Bulge_Size      Total   PAM_gen Var_uniq        Samples Annotation Type Real Guide
#NOTE1 can happend that a iupac falls under the N char of NGG, meaning that a target can have the same number of mms both in his REF and VAR part:
#CACTGCAACCTCTGTCTCCCKGG
#CACTGCAACCTCTGTCTCCCGGG        REF
#CACTGCAACCTCTGTCTCCCTGG        VAR
#So check if decrease_ref_count is not empty to avoid this (in this case +1 will be added to samples for VAR part of target and +1 for all the
# other samples for REF part)

#Block scikit warnings
import warnings 
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings('ignore',category=UserWarning)

import sys
import json
import time
import itertools
import os
from os.path import isfile, isdir,join      #for getting lst of chr to know file extension and if enriched
from os import listdir
from intervaltree import Interval, IntervalTree
import concurrent.futures
import subprocess
import pandas as pd
import pickle       #to read CFD matrices
import numpy as np
import azimuth.model_comparison
import string
import multiprocessing
import re
from supportFunctions.loadSample import associateSample


SIZE_DOENCH = 10000
N_THR = 3
step = sys.argv[-1] #IF 'Step' in step then do prints to show the user the current analysis status

#Return max doench value among list of extended targets
def doenchParallel(targets, model, result):
    doench_score =  azimuth.model_comparison.predict(targets,None, None, model= model, pam_audit=False)
    doench_score = [np.around(i * 100) for i in doench_score]
    max_doench = int(max(doench_score))
    result.append(max_doench)

def doenchForIupac(sequence_doench, guide_seq, genome_type): 
    pos_iupac = []
    var = []
    for pos, c in enumerate(sequence_doench):
        if c in iupac_code:
            pos_iupac.append(pos)
            var.append(iupac_code[c])
  
    if var:
        for i in itertools.product(*var):
            t = list(sequence_doench)
            for p, el in enumerate(pos_iupac):
                t[el] = i[p]
            targets_for_doench[guide_seq][genome_type].append(''.join(t))
    else:
        targets_for_doench[guide_seq][genome_type].append(sequence_doench)

def get_mm_pam_scores():
    try:
        mm_scores = pickle.load(open(os.path.dirname(os.path.realpath(__file__)) + '/mismatch_score.pkl', 'rb'))
        pam_scores = pickle.load(open(os.path.dirname(os.path.realpath(__file__)) +'/PAM_scores.pkl', 'rb'))
        return (mm_scores, pam_scores)
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")


def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', '-':'-'}
    letters = list(s[::-1])
    try:
        letters = [basecomp[base] for base in letters]
    except:
        return None     #If some IUPAC were not translated
    return ''.join(letters)

# Calculates CFD score
def calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores):
    score = 1
    sg = sg.replace('T', 'U')
    guide_seq = guide_seq.replace('T', 'U')
    s_list = list(sg)
    guide_seq_list = list(guide_seq)
    for i, sl in enumerate(s_list):  
        if guide_seq_list[i] == sl:
            score *= 1
        else:
            try:        #Catch exception if IUPAC characterv
                key = 'r' + guide_seq_list[i] + ':d' + revcom(sl) + ',' + str(i + 1)
            except:
                score = 0
                break
            try:
                score *= mm_scores[key]
            except: #If '-' is in first position, i do not have the score for that position
                pass
    score *= pam_scores[pam]
    return score
if 'Step' in step:
    print(step + ': Reading Input Files', end = '\r')
else:
    print("READING INPUT FILES")
#Dictionaries for annotating samples

#Dict for populations
# pop_file = pd.read_excel(os.path.dirname(os.path.realpath(__file__)) + '/20130606_sample_info.xlsx')
# all_samples = pop_file.Sample.to_list()
# all_pop = pop_file.Population.to_list()
# dict_sample_to_pop = dict()
# for  pos, i in enumerate(all_samples):
#     try:
#         dict_sample_to_pop[i] = all_pop[pos]        #{'S1':'POP1', 'S2':'POP1', ...}
#     except:
#         dict_sample_to_pop[i] = all_pop[pos]

# #Dict for superpopulation
# dict_pop_to_sup = {'CHB':'EAS', 'JPT':'EAS', 'CHS':'EAS', 'CDX':'EAS', 'KHV':'EAS',
#                     'CEU':'EUR', 'TSI':'EUR', 'FIN':'EUR', 'GBR':'EUR', 'IBS':'EUR',
#                     'YRI':'AFR', 'LWK':'AFR', 'GWD':'AFR', 'MSL':'AFR', 'ESN':'AFR', 'ASW':'AFR', 'ACB':'AFR',
#                     'MXL':'AMR', 'PUR':'AMR', 'CLM':'AMR', 'PEL':'AMR',
#                     'GIH':'SAS', 'PJL':'SAS', 'BEB':'SAS', 'STU':'SAS', 'ITU':'SAS'
# }
# superpopulation = ['EAS', 'EUR', 'AFR', 'AMR','SAS']

dict_sample_to_pop, dict_pop_to_sup, dict_superpop_to_pop, dict_pop_to_sample, all_samples, all_pop, superpopulation, gender_sample = associateSample.loadSampleAssociation(sys.argv[11])


#READ INPUT FILES
annotationFile = sys.argv[1] #file with annotation
resultsFile = sys.argv[2] #file with results from search
outputFile = sys.argv[3] #file with annotated results

#Get pam and guide length for new count mismatch samples
pam_at_beginning = False
with open (sys.argv[5]) as pam:
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

outFileCFD = open(outputFile + '.targets.CFD.txt','w')# file with targets and the corresponding CFD score
do_scores = True
if guide_len != 20 or len(pam) !=3 or pam_at_beginning:
    with open(outputFile + '.scores.txt', 'w+') as result:
        result.write('NO SCORES')
        outFileCFD.write('NO SCORES')
        do_scores = False

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


#OPEN INPUT FILES AND PREPARE OUTPUT FILE
inResult = open(resultsFile, "r")  # resultfile open
inAnnotationFile = open(annotationFile, "r")  # file with annotations open
# outFileSampleAll = open(outputFile + '.samples.all.annotation.txt', 'w')  # outfile open (file with IUPAC targets and associated samples and annotation)
outFileSample = open(outputFile + '.samples.annotation.txt', 'w') #file with real nucleotides with associated samples and annotation
outFileSummary = open(outputFile + '.Annotation.summary.txt', 'w')  # outfile open (summary file calculated on top1file)

process = subprocess.Popen(['wc', '-l', resultsFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = process.communicate()
total_line = int(out.decode('UTF-8').split(' ')[0])
if total_line < 1:
    print('WARNING! Input file has no targets')
    sys.exit()
if total_line < 10:
    mod_tot_line = 1
else:
    mod_tot_line = int(total_line/10)
#VARIABLE INIT
guideDict = {}
totalDict = {}

start_time = time.time()

if 'Step' in step:
    print(step + ': Executing Preliminary Operations', end = '\r')
else:
    print("EXECUTING PRELIMINARY OPERATIONS")

annotationsTree = IntervalTree()
annotationsSet = set()
#guidesSet = set()       #NOTE/BUG if guide finds 0 targets, it will not be annotated

for line in inAnnotationFile:
    x = line.split('\t')
    x[3] = str(x[3]).rstrip("\n")
    annotationsTree[int(x[1]):int(x[2])] = str(x[0])+'\t'+str(x[3])
    annotationsSet.add(str(x[3]))

totalDict['targets'] = [0]*10
for item in annotationsSet:
    totalDict[item] = [0]*10

if 'Step' not in step:
    print("PRELIMINARY OPERATIONS COMPLETED IN: %s seconds" % (time.time() - start_time))

start_time = time.time()
if 'Step' in step:
    print('\033[K' + step + ': Executing Samples Extraction and Annotation', end = '\r')
else:
    print("EXECUTING SAMPLE PROCESSING AND ANNOTATION")

with open(resultsFile, 'r') as resFile:
    header_len = len(resFile.readline().strip().split('\t'))

# if header_len == 15:    #'Both' case : comparison variant/ref is active
header = '#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation Type\tReal Guide'


mm_pos = 7      #position of mismatch column
bulge_pos = 8
max_dna_bulges = int(sys.argv[9])
max_rna_bulges = int(sys.argv[10])
max_bulges = max_dna_bulges
if max_rna_bulges > max_bulges:
    max_bulges = max_rna_bulges

blank_add_begin = ' '               #Needed for replacing IUPAC in cluster targets
blank_add_end = ''
pam_multiplier = 1
pam_multiplier_negative = 0
start_sample_for_cluster = 0
cluster_step = 1    #If PAM end, go left to right
sum_for_mms = 0 #when updatig lowercase for nem_mm, this value represents the offset for the pam position (mainly needed only if pam at beginning)
end_sample_for_cluster = max_dna_bulges + max_rna_bulges  #Values to check new iupac when working on cluster targets
if pam_at_beginning:
    blank_add_begin = ''
    blank_add_end = ' '
    pam_multiplier = 0              #Since ' ' are at end, and '-' to reinsert are before the ' ', need to put max_dna_bulges and rna_bulges of target to 0
    pam_multiplier_negative = 1
    end_sample_for_cluster = len_pam + guide_len - max_rna_bulges
    start_sample_for_cluster = len_pam + guide_len + max_dna_bulges
    cluster_step = -1   #If PAM beginning, go right to left
    sum_for_mms = len_pam
outFileSample.write(header + '\n')
if do_scores:
    outFileCFD.write(header + '\tCFD\n')
# outFileSampleAll.write(header + '\n')
summary_samples = True

header_list = header.strip().split('\t')
#Variables for summary samples code
'''
{
    GUIDE1 -> {
        SAMPLE/POP/SUPERPOP1 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        },
        SAMPLE/POP/SUPERPOP2 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        }
    }
    GUIDE2 -> {
        SAMPLE/POP/SUPERPOP1 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        },
        SAMPLE/POP/SUPERPOP2 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        }
    }
}

Per pop e superpop, se ho due sample stessa famiglia stesso target, conto solo una volta (visited_pop and visited_superpop array)
'''
count_sample = dict()       #NOTE cout_sample -> GUIDE -> SAMPLE -> has targets + ann1 + ann2 ... + refposition
    #refposition is a key unrelated to the other keys (targets ann1 ann2 ...) and it's used to classify the sample (0 0+ 1 1+).
    #it's put in here just to avoid to duplicate the entire guide -> sample ->     structure
    # refposition -> [class , number of specific VAR on target to add/remove]  #Save class (0 at start) and number of ontarget var specific
    #for that sample.
ontarget_reference_count = dict() #Count number of REF target and REF part of semicommon
count_pop = dict()
count_superpop = dict()     #NOTE added key 'distributions' for population distribution images
    #count_superpop-> GUIDE -> SUPERPOP -> targets ann1 ann2 ... distributions
    # distributions is an array of len mms+bulge, each position contains an array [0,0,0] of len bulge+1 (indicating no bulge, 1 bulge, 2bulge ...)

#Create -Summary_total for a file ref.Annotation.summary.txt from the y and n values of Var_uniq column
summary_barplot_from_total = False
if 'Var_uniq' in header:
    summary_barplot_from_total = True
    vu_pos = header_list.index('Var_uniq')
count_unique = dict()
count_unique['targets'] = [0]*10
count_unique_for_guide = dict()
for item in annotationsSet:
    count_unique[item] = [0]*10

#Variables for samples calculation
total_error = 0


current_chr = 'none'
chr_name = 'none'

def rev_comp(a):
    if a == 'A' or a == 'a':
        return 'T'
    if a == 'T' or a == 't':
        return 'A'
    if a == 'C' or a == 'c':
        return 'G'
    return 'C'

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

#For scoring of CFD And Doench
tab = str.maketrans("ACTGRYSWMKHDBVactgryswmkhdbv", "TGACYRSWKMDHVBtgacyrswkmdhvb") 

def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]

mm_scores, pam_scores = get_mm_pam_scores()
guides_dict = dict()
guides_dict_doench = dict()
targets_for_doench = dict()

N_THR = multiprocessing.cpu_count() // 2
refgenomedir = sys.argv[7]
add_enr = ''
list_chr = [f for f in listdir(refgenomedir) if isfile(join(refgenomedir, f)) and not f.endswith('.fai')] 
if '.enriched.' in list_chr[0]:
    add_enr = '.enriched'
elif '.indels.' in list_chr[0]:
    add_enr = '.indels'
add_ext = '.fa'
if '.fasta' in list_chr[0]:
    add_ext = '.fasta'

with open( os.path.dirname(os.path.realpath(__file__)) + "/azimuth/saved_models/V3_model_nopos.pickle", 'rb') as f:
    model = pickle.load(f)
max_doench = 0
sum_cfd = 0
cfd_scores = []

remove_tmp_bedfile = False
start_time_total = time.time()
lines_processed = 0
allowed_mms = int(sys.argv[6])
current_guide_chr_pos = 'no'
cluster_update = open(outputFile + '.cluster.tmp.txt', 'w+')
cluster_update.write(header + '\n')        #Write header
save_cluster_targets = True
remove_iupac = False
save_total_general_table = False
add_to_general_table = dict()   #chiave ['add'] = target semicommon da aggiungere alla tab generale delle guide, metto i valori di total presi dal primo target REF nel cluster di un semicommon che esiste
                                #chiave ['distribution'] = array len total, ogni cella divisa per mm,1B,2B..., per populationDistribution
last_annotation = ''    #needed when counting the ref part of a semicommon in order to not redo the annotation
last_samples = ''       #needed for save cluster part
next(inResult)      #Skip header
for line in inResult:
    x = line.strip().split('\t')
    guide_no_bulge = x[1].replace("-","")
    if (guide_no_bulge + x[3] + x[5] + x[6]) == current_guide_chr_pos:     #Target is in current cluster, simply save the sample and annotation, discard if status is F
        if save_cluster_targets:
            if remove_iupac:        #Save only Reference target of the cluster
                for c in x[2]:
                    if c in iupac_code:
                        break
                else:       #no break triggered
                    x[13] = last_annotation 
                    cluster_update.write('\t'.join(x) + '\n')       #Save the REF target, update only the Annotation
                continue        
            #Keep the semicommon ref total value to be added to general table
            if save_total_general_table:
                for c in x[2]:
                    if c in iupac_code:
                        break
                else: #Found first REF target in cluster of semicommon
                    if x[bulge_pos + 1] == '0':     #If total column is 0  -> On-target
                        ontarget_reference_count[guide_no_bulge] +=1        #Add +1 to Ref count
                        if decrease_ref_count: #VIEW NOTE1
                            for s in decrease_ref_count.split(','):                        #But Decrease count for samples that do not have the REF (ie sample with ok target but not 0 total, or sample with not ok target)
                                try:
                                    count_sample[guide_no_bulge][s]['refposition'][2] -= 1
                                except:
                                    count_sample[guide_no_bulge][s]['refposition'] = [0,0,-1]
                                count_sample[guide_no_bulge][s]['refposition'][0] += 1      #Visited +1
                         
                    x[13] = last_annotation
                    outFileSample.write('\t'.join(x) + '\n')      
                    add_to_general_table[guide_no_bulge]['add'][int(x[mm_pos]) +  int(x[bulge_pos])] += 1  
                    #Do annotation to keep numbers consistent between images and general table
                    #conto i target generali per mm threshold
                    totalDict['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                    guideDict[guide_no_bulge]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                    #REF target, add +1 to PopulationDistributionREF
                    add_to_general_table[guide_no_bulge]['distributions'][int(x[mm_pos]) + int(x[bulge_pos])][int(x[bulge_pos])] +=1   #add +1 to array in cell TOTAL, in comma position corresponding to bulge size
                    
                    #Conto per annotazione
                    for ann in last_annotation.split(','):
                        if ann == 'n':
                            break
                        guideDict[guide_no_bulge][ann][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                        totalDict[ann][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                    
                    #Calculate scores
                    if do_scores:       #Calculate scores for reference targets
                        if x[0] == 'DNA':
                            cfd_score = calc_cfd(x[1][int(x[bulge_pos]):], x[2].upper()[int(x[bulge_pos]):-3], x[2].upper()[-2:], mm_scores, pam_scores)
                        else:
                            cfd_score = calc_cfd(x[1], x[2].upper()[:-3], x[2].upper()[-2:], mm_scores, pam_scores)
                        outFileCFD.write('\t'.join(x) + '\t' + str(cfd_score) + '\n')
                        sum_cfd = sum_cfd + cfd_score
                        try:
                            guides_dict[x[1]] = guides_dict[x[1]] + cfd_score
                        except:
                            guides_dict[x[1]] = cfd_score

                        if x[bulge_pos + 1] == '0':    #DOENCH
                            #estraggo sequenza
                            with open(outputFile + '.bedfile_tmp.bed', 'w+') as bedfile:
                                remove_tmp_bedfile = True
                                if x[6] == '+':
                                    bedfile.write(x[3] + '\t' + str(int(x[4]) - 4 ) + '\t' + str(int(x[4]) + 23 + 3 ))
                                else:
                                    bedfile.write(x[3] + '\t' + str(int(x[4]) - 3 ) + '\t' + str(int(x[4]) + 23 + 4 ))
                            #Extract sequence from REFERENCE
                            extr = subprocess.Popen(['bedtools getfasta -fi ' + refgenomedir + '/' + x[3] + add_enr + add_ext + ' -bed ' + outputFile + '.bedfile_tmp.bed'], shell = True, stdout=subprocess.PIPE)  #TODO insert option for .fasta
                            extr.wait()
                            out, err = extr.communicate()
                            out = out.decode('UTF-8')
                            if x[6] == '+':
                                sequence_doench = out.strip().split('\n')[-1].upper()
                                # sequence_doench = sequence_doench[:4] + x[2] + sequence_doench[-3:]
                            else:
                                sequence_doench = reverse_complement_table(out.strip().split('\n')[-1].upper())
                                # sequence_doench = sequence_doench[:4] + x[2] + sequence_doench[-3:]
                            
                            if x[1] not in targets_for_doench:
                                targets_for_doench[x[1]] = {'ref': [], 'enr': []}
                            doenchForIupac(sequence_doench, x[1], 'ref')  #Get all possible targets with iupac itertools for doench
                    save_total_general_table = False
            #Save target (with or without IUPAC in cluster file)
            #Put top1 scomposed nucleotide in target in cluster that have IUPAC characters
            line = line.strip().split('\t')
            for c in line[2]:
                if c in iupac_code:
                    break
            else:   #No break triggered: target doesn't have iupac
                line[13] = last_annotation
                cluster_update.write('\t'.join(line) + '\n') 
                continue
            #From now, all the target have at least one IUPAC
            tmp_gap_position = []
            if line[0] == 'X':  
                target_to_modify = list(blank_add_begin * max_dna_bulges + line[2] + blank_add_end * max_dna_bulges)        #Eg if max bulge dna is 1: pam end -> [' ', 'A','C','G',...,'A','G','G']; pam begin -> ['T','G','T','C','A','G',...,'A','T',' '] 
            elif line[0] == 'DNA':
                target_to_modify = list(blank_add_begin * (max_dna_bulges - int(line[bulge_pos])) + line[2] + blank_add_end * (max_dna_bulges - int(line[bulge_pos])))
            else:
                tmp_gap_position = [g.start() for g in re.finditer('-', line[2])]
                target_to_modify = list(blank_add_begin * (max_dna_bulges + int(line[bulge_pos])) + line[2].replace('-','') + blank_add_end * (max_dna_bulges + int(line[bulge_pos])))   
            
            for elem in  pos_char_for_cluster:      #elem is tuple (position in target, character that was substituted in top1scomposed)
                if target_to_modify[elem[0]] == ' ':       #if target to modify is ' '
                    continue
                if target_to_modify[elem[0]] not in iupac_code:                 #Non IUPAC char should not be targeted
                    print('Error: Substituting into a NON IUPAC character:', line, str(elem[0]) , pos_char_for_cluster) 
                target_to_modify[elem[0]] = elem[1]
                
            #Fix the removed '-' in RNA targets
            for tmp_g_p in tmp_gap_position:
                target_to_modify.insert(tmp_g_p + (max_dna_bulges * pam_multiplier) + int(line[bulge_pos]) * pam_multiplier, '-')
            
            #DNA and X cases can have more characters at beginning
            remove_target = False
            iupac_pos = 0
            bulge_found = 0
            target_to_modify = [c for c in target_to_modify if c != ' ' ]
            total_bulges = 0
            pam_sum_cluster = 0
            if line[0] == 'RNA':
                total_bulges = int(line[bulge_pos])
                pam_sum_cluster = 1
            tmp_intersection = set()
            for i in range (min(start_sample_for_cluster, len(target_to_modify) - 1), end_sample_for_cluster, cluster_step):    #Check only on some characters, not all target
                if target_to_modify[i] == '-':
                    bulge_found += 1
                if target_to_modify[i] in iupac_code:
                    if line[6] == '-':
                        # iupac_pos = str(int(line[4]) + (len(target_to_modify) - i - 1) + 1 - (no_char * pam_multiplier_negative)  )
                        iupac_pos = str(int(line[4]) + (len(target_to_modify) - i) - (total_bulges - bulge_found) + pam_sum_cluster * pam_multiplier_negative)
                    else:
                        # iupac_pos = str(int(line[4]) + i + 1 - no_char)
                        iupac_pos = str(int(line[4]) + i + 1 - (total_bulges - bulge_found) * pam_multiplier_negative)
                    try:
                        a = (datastore[chr_name + ',' + iupac_pos])   #NOTE se non ha samples, ritorna ;ref,var
                    except Exception as e:
                        print(e)
                        print('Error at', target_to_modify, ' pos line', line[4], ', iupac pos', str(iupac_pos), ', i:', str(i))
                        continue
                    else:
                        a = a.split('/')
                        for samples_chars in a:     #samples_char can be 'sample1,sample2;A,T' or ';A,T'
                            samples_chars = samples_chars.split(';')
                            ref_char = samples_chars[-1].split(',')[0]
                            var_char = samples_chars[-1].split(',')[1]
                
                            if line[6] == '-':
                                ref_char = rev_comp(ref_char)
                                var_char = rev_comp(var_char)

                            if 'n' in last_samples:
                                target_to_modify[i] = ref_char  #RefChar
                                break       #Go to next iupac

                            #Some samples were already calculated for the target
                            sample_set = set(samples_chars[0].split(','))   #set of sample associated with the current IUPAC
                            tmp_intersection = sample_set & last_samples        #Intersection
                            if tmp_intersection:
                                if len(last_samples) == len(tmp_intersection):      #Put var char only if all sample of the top1 are in this new iupac
                                    target_to_modify[i] = var_char   #VarChar
                                    break   #Go to next iupac
                                else:
                                    remove_target = True
                                    break   #There's at least one target in common, but in this else not all samples are in common, meaning  that the next sample_chars either has no common sample or <len(lastsamples) since it will not have this common sample
                        else: #No break encountered #No intersection on all the samples_chars, put ref char
                            target_to_modify[i] = ref_char   #RefChar

            if remove_target:
                continue
            else:       #Do also a check for new mm value and if mm value > threshold, do not save target
                line[2] = ''.join(target_to_modify).strip()
                mm_new_t = 0
                guide_no_pam = line[1][pos_beg:pos_end]    
                list_t = list(line[2])
                for position_t, char_t in enumerate(line[2][pos_beg:pos_end]): 
                    if char_t.upper() != guide_no_pam[position_t]:
                        mm_new_t += 1
                        if guide_no_pam[position_t] != '-':
                            list_t[sum_for_mms + position_t] = char_t.lower()
                    
                if allowed_mms < (mm_new_t - int(line[bulge_pos])):        
                    continue                #Remove target since id does not respect mms constrains
                line[2] = ''.join(list_t)
                line[mm_pos] = str(mm_new_t - int(line[bulge_pos]))
                line[bulge_pos + 1] = str(mm_new_t) #total differences between targets and guide (mismatches + bulges)
            
            line[12] = ','.join(last_samples)
            line[13] = last_annotation
            line = '\t'.join(line) 
            cluster_update.write(line + '\n') 
        # lines_processed +=1
        # if lines_processed % (mod_tot_line) == 0:
        #     print('Annotation: Total progress ' + str(round(lines_processed /total_line *100, 2)) + '%')
        continue
    save_cluster_targets = True
    remove_iupac = False
    decrease_ref_count = []     #Samples that does not have REF or VAR part of semicommon
    pos_char_for_cluster = []       #contains (pos, char) where to change nucleotides for writing the cluster targets --> see conversione_caratteri_cluster.txt
    current_guide_chr_pos = guide_no_bulge + x[3] + x[5] + x[6]
    if x[3] != current_chr:
        if not os.path.exists(os.path.realpath(sys.argv[4]) + '/my_dict_' + x[3] + '.json'):
            pass
        else:
            # if 'none' not in current_chr:
            #     print('Done ', current_chr)
            current_chr = x[3]
            chr_name = x[3]
            with open(os.path.realpath(sys.argv[4]) + '/my_dict_' + current_chr + '.json', 'r') as f:
                start_time = time.time()
                datastore = json.load(f)
                if 'Step' in step:
                    print('\033[K' + step + ': Processing ' + current_chr, end = '\r')
                else:
                    print ('PROCESSING ' + current_chr)
    
    pos_snp = []
    tuple_var_ref = []
    target_combination = []
    pos_snp_chr = []
    set_list = []
    target_string = x[2]
    if x[6] == '-':
        target_string = target_string[::-1]
    bulge_found = 0 
    for pos, char in enumerate(target_string):
        if char == '-':
            bulge_found = bulge_found + 1 
        if char in iupac_code:
            iupac_pos = str(int(x[4]) + pos + 1 - bulge_found)
            try:
                a = (datastore[chr_name + ',' + iupac_pos])   #NOTE se non ha samples, ritorna ;ref,var
            except Exception as e:      #NOTE this error can occure if i have an IUPAC in a target that has no vcf file
                print(e)
                print('Error at ' + line.rstrip() + ', with char ' + char + ', at pos ', iupac_pos, '. No corresponding SNP position was found in the vcf file')
                samples_this_position = []
                total_error = total_error + 1
            else:   
                a = a.split('/')
                samples_this_position = []      #Set of samples in this position (even if they have different var character)
                character_list =  []          #List of var characters, later add the reference character
                for samples_chars in a:     #samples_char can be 'sample1,sample2;A,T' or ';A,T'
                    samples_chars = samples_chars.split(';')
                    ref_char = samples_chars[-1].split(',')[0]
                    var_char = samples_chars[-1].split(',')[1]
                    if x[6] == '-':
                        ref_char = rev_comp(ref_char)
                        var_char = rev_comp(var_char)
                    character_list.append(var_char)
                    if samples_chars[0]:
                        samples_this_position.extend(samples_chars[0].split(','))
                character_list.append(ref_char)    
                pos_snp.append(pos)
                pos_snp_chr.append(iupac_pos)
                tuple_var_ref.append(tuple(character_list))
            if samples_this_position:
                set_list.append(set(samples_this_position))
            else:
                set_list.append(set())
    #Get Union of all samples
    union_sample = list(set().union(*set_list))
    if union_sample:
        x[12] = ','.join(union_sample)
    else:
        x[12] = 'n'

    #Get all combinations to remove targets that have no haplotype 
    #Create all combinations
    for i in itertools.product(*tuple_var_ref):
        t = list(target_string)
        for p, el in enumerate(pos_snp):
            t[el] = i[p]
        target_combination.append(''.join(t))
    
    target_scomposti_salvare = []
    samples_already_assigned = set()
    false_targets = 0

    #Get PAM of reference
    if x[6] == '-':
        tmp_t = target_combination[-1][::-1]
        reference_pam = tmp_t[pam_begin:pam_end]
    else:
        reference_pam = target_combination[-1][pam_begin:pam_end]
    
    found_creation = False
    for pos_pam_ref, pam_char_ref in enumerate(reference_pam):
        if not iupac_code_set[pam[pos_pam_ref]] & iupac_code_set[pam_char_ref]:     #ref char not in set of general pam char
            found_creation = True

    for t in target_combination:
        set_list2 = []
        final_result = x.copy()
        for ele_pos,p in enumerate(pos_snp_chr):
            try:
                a = (datastore[chr_name + ',' + p])        
            except:
                set_list2.append(set())
            else:
                a = a.split('/')
                for samples_chars in a:     #samples_char can be 'sample1,sample2;A,T' or ';A,T'
                    samples_chars = samples_chars.split(';')
                    # ref_char = samples_chars[-1].split(',')[0]
                    var_char = samples_chars[-1].split(',')[1]
                    if x[6] == '-':
                        # ref_char = rev_comp(ref_char)
                        var_char = rev_comp(var_char)
                    
                    if t[pos_snp[ele_pos]].upper() == var_char:
                        if samples_chars[0]:
                            set_list2.append(set(samples_chars[0].split(',')))
                        else:
                            set_list2.append(set())
                        break       
            
        if set_list2:
            common_samples = set.intersection(*set_list2)
            common_samples = common_samples - samples_already_assigned
            samples_already_assigned = samples_already_assigned.union(common_samples)
            if common_samples:
                final_result[12] = ','.join(common_samples)
            else:
                # final_result.append('No common samples')
                final_result = []                       #DO not save results without samples
                false_targets += 1
        else:
            # final_result.append('No samples')         #DO not save results without samples
            final_result = []
            if set_list:            #Increase false_targets on targets that have at least 1 IUPAC
                false_targets += 1
        if x[6] == '-':
            t = t[::-1]
        mm_new_t = 0
        
        if final_result:
            guide_no_pam = final_result[1][pos_beg:pos_end]    
            list_t = list(t)  
            for position_t, char_t in enumerate(t[pos_beg:pos_end]):
                if char_t.upper() != guide_no_pam[position_t]:
                    mm_new_t += 1
                    if guide_no_pam[position_t] != '-':
                        list_t[sum_for_mms + position_t] = char_t.lower()
            final_result[2] = ''.join(list_t)#t
        
            #Check for pam status
            pam_ok = True
            for pam_chr_pos, pam_chr in enumerate(t[pam_begin:pam_end]):
                if pam_chr.upper() not in iupac_code_set[pam[pam_chr_pos]]:
                    pam_ok = False

            if not pam_ok or allowed_mms < (mm_new_t - int(final_result[8])):                     
                false_targets += 1
                x[12] = ','.join(set(x[12].split(',')) - set(common_samples))
                decrease_ref_count.append(','.join(common_samples))   #This samples does not have the VAR version nor the REF version
                continue                #Remove target since id does not respect PAM or mms constrains

            final_result[7] = str(mm_new_t - int(final_result[8]))
            final_result[9] = str(mm_new_t) #total differences between targets and guide (mismatches + bulges)
            if found_creation:
                final_result[10] = t[pam_begin:pam_end]
            
            target_scomposti_salvare.append(final_result)
    if false_targets >= len(target_combination):        #If all the scomposed targets have no sample or do not respect PAM/mm threasold, the iupac target does not really exist
        line = line.strip().split('\t')
        if line[vu_pos] == 'y':         #Do not do annotation because target does not exists, and do not save his cluster
            save_cluster_targets = False
            continue    #DO NOT save this target because no ref homologous and no sample combination exists
        
        x[12] = 'n'     #Since iupac target has no scomposition, it means it has no sample associated
        line = '\t'.join(line)
        save_cluster_targets = True        #salvare il cluster ma togliendo gli iupac
        remove_iupac = True
        x = next(inResult).strip().split('\t')   #get next target of the cluster
        entered = False
        while (x[1].replace('-','') + x[3] + x[5] + x[6]) == current_guide_chr_pos:   #while still in same cluster
            entered = True
            for c in x[2]:
                if c in iupac_code:
                    break
            else:   #no break triggered in previous for --> x[2] has no iupac char
                break
            x = next(inResult).strip().split('\t')
        if not entered:
            continue 
        line = '\t'.join(x)                     #Fist target in the cluster that is REF
        # x.append('n')                       #Fist target in the cluster that is REF -> add 'n' as sample
        # x.append(x[1].replace('-',''))      #Fist target in the cluster that is REF -> add real guide
        tuple_var_ref = []      #Since this is now a REF target, it has no iupac --> needed to save to sample.annotation file
    
    
    if target_scomposti_salvare:        #Keep the target with lowest total and mms as representative of the IUPAC target
        if x[vu_pos] == 'n':                #Save mms of REF of a semicommon
            save_total_general_table = True
        target_scomposti_salvare.sort(key = lambda x : (int(x[mm_pos + 2]), int(x[mm_pos]))) #Order scomposition by total and mms values
        x[2] = target_scomposti_salvare[0][2]       #Adjust Target sequence, from IUPAC to first of scomposition
        x[mm_pos] = target_scomposti_salvare[0][mm_pos]
        x[mm_pos + 2] = target_scomposti_salvare[0][mm_pos + 2]     #Adjust IUPAC with min total and mms of his scomposition
        original_targ_len = len(target_scomposti_salvare[0][2])
        if target_scomposti_salvare[0][0] == 'RNA':
            gap_position = [g.start() for g in re.finditer('-', target_scomposti_salvare[0][2])]    #list of indices where '-' is located
        for elem in pos_snp:
            add_to_count = 0
            add_blank = int(target_scomposti_salvare[0][bulge_pos])
            if target_scomposti_salvare[0][0] == 'RNA':
                add_blank = 0
                for i in gap_position:
                    if not pam_at_beginning:        #TODO try to remove the pam_at_beginning check to improve performances (it's done only on top1 so maybe it's ok)
                        if target_scomposti_salvare[0][6]=='-':
                            if (original_targ_len - elem - 1) < i:
                                add_to_count += 1
                        else:
                            if elem < i:
                                add_to_count += 1
                    else:
                        if target_scomposti_salvare[0][6]=='-':
                            if (original_targ_len - elem - 1) > i:
                                add_to_count -= 1
                        else:
                            if elem > i:
                                add_to_count -= 1

            if target_scomposti_salvare[0][6]=='-':
                pos_char_for_cluster.append(((len(target_scomposti_salvare[0][2]) - elem - 1) + add_to_count + (max_dna_bulges - int(add_blank) ) * pam_multiplier, target_scomposti_salvare[0][2][(len(target_scomposti_salvare[0][2]) - elem - 1)]))
            else:
                pos_char_for_cluster.append((elem + add_to_count + (max_dna_bulges - int(add_blank) ) * pam_multiplier, target_scomposti_salvare[0][2][elem]))   #save (position, character) of best scomposition
    
    #Annotate target
    visited_pop = []
    visited_superpop = []

    #inserisco la key nel dict se non presente e creo la sua matrice
    if(guide_no_bulge not in guideDict.keys()):
        guideDict[guide_no_bulge] = {}
        guideDict[guide_no_bulge]['targets'] = {}
        guideDict[guide_no_bulge]['targets'] = [0]*10

        add_to_general_table[guide_no_bulge] = dict()
        add_to_general_table[guide_no_bulge]['add'] = [0] * 10    # GUIDE -> ['add'] -> [ 0 0 0 0 0 ...] valori per total (mms + bulge)
        add_to_general_table[guide_no_bulge]['distributions'] = [[0] * (max_bulges + 1) for x in range(10)] # GUIDE -> ['distributions'] -> [ 0,0,0 0,0,0 0,0,0 ...] len array = 10, ogni cella ha tanti 0 quanti numero di max_bulge +1

        count_unique_for_guide[guide_no_bulge] = dict()                 #NOTE count_unique means that the target have at least 1 sample
        count_unique_for_guide[guide_no_bulge]['targets'] = [0]*10

        count_sample[guide_no_bulge] = dict()
        count_pop[guide_no_bulge] = dict()
        count_superpop[guide_no_bulge] = dict()
        ontarget_reference_count[guide_no_bulge] = 0 #List of chr,clusterpos of on target of reference, needed for classify samples
        for item in annotationsSet:
            guideDict[guide_no_bulge][item]= {}
            guideDict[guide_no_bulge][item] = [0]*10

            count_unique_for_guide[guide_no_bulge][item] = [0]*10
    
    #conto i target generali per mm threshold
    totalDict['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
    guideDict[guide_no_bulge]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1

    if summary_barplot_from_total:
        if x[12] != 'n':            #If the target has at least 1 sample, initialize count_unique
            count_unique['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            count_unique_for_guide[guide_no_bulge]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
    
    if summary_samples:
        for sample in x[12].split(','):
            if sample == 'n':
                #REF target, add +1 to PopulationDistributionREF
                add_to_general_table[guide_no_bulge]['distributions'][int(x[mm_pos]) + int(x[bulge_pos])][int(x[bulge_pos])] +=1   #add +1 to array in cell TOTAL, in comma position corresponding to bulge size
                continue
            #Initialization if sample, pop or superpop not in dict
            if sample not in count_sample[guide_no_bulge]:
                count_sample[guide_no_bulge][sample] = {'targets': [0]*10}
                for item in annotationsSet:
                    count_sample[guide_no_bulge][sample][item] = [0]*10
                count_sample[guide_no_bulge][sample]['refposition'] = [0, 0, 0]        #[Visited count, new ontargetVAR, noOntargetVAR/REF]Needed for classify samples
                if dict_sample_to_pop[sample] not in count_pop[guide_no_bulge]:
                    count_pop[guide_no_bulge][dict_sample_to_pop[sample]] = {'targets': [0]*10}
                    for item in annotationsSet:
                        count_pop[guide_no_bulge][dict_sample_to_pop[sample]][item] = [0]*10
                if dict_pop_to_sup[dict_sample_to_pop[sample]] not in count_superpop[guide_no_bulge]:
                    count_superpop[guide_no_bulge][dict_pop_to_sup[dict_sample_to_pop[sample]]] = {'targets': [0]*10}
                    for item in annotationsSet:
                        count_superpop[guide_no_bulge][dict_pop_to_sup[dict_sample_to_pop[sample]]][item] = [0]*10
                    count_superpop[guide_no_bulge][dict_pop_to_sup[dict_sample_to_pop[sample]]]['distributions'] = [[0] * (max_bulges + 1) for x in range(10)]
            #Add +1 to targets
            count_sample[guide_no_bulge][sample]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            if dict_sample_to_pop[sample] not in visited_pop:
                visited_pop.append(dict_sample_to_pop[sample])
                count_pop[guide_no_bulge][dict_sample_to_pop[sample]]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            if dict_pop_to_sup[dict_sample_to_pop[sample]] not in visited_superpop:
                visited_superpop.append(dict_pop_to_sup[dict_sample_to_pop[sample]])
                count_superpop[guide_no_bulge][dict_pop_to_sup[dict_sample_to_pop[sample]]]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                count_superpop[guide_no_bulge][dict_pop_to_sup[dict_sample_to_pop[sample]]]['distributions'][int(x[mm_pos]) + int(x[bulge_pos])][int(x[bulge_pos])] +=1 #add +1 to array in pos TOTAL, in cell corrresponding to bulge size
        visited_pop = []
        visited_superpop = []
    
    #faccio match su albero
    foundAnnotations = sorted(annotationsTree[int(x[4]):(int(x[4])+int(len(guide_no_bulge))+1)])
    string_annotation = []
    found_bool = False
    for found in range(0, len(foundAnnotations)):
        guide = foundAnnotations[found].data
        guideSplit = guide.split('\t')
        # print(guide, str(guideSplit[0]), str(x[3]))
        if(str(guideSplit[0]) == str(x[3])):
            found_bool = True
            #outFileTargets.write(line.rstrip() + '\t' + str(guideSplit[1]) + "\n")
            string_annotation.append(str(guideSplit[1]))
            guideDict[guide_no_bulge][guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            totalDict[guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            
            if summary_barplot_from_total:
                if x[12] != 'n':
                    count_unique[guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                    count_unique_for_guide[guide_no_bulge][guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            
            if summary_samples:
                for sample in x[12].split(','):
                    if sample == 'n':
                        continue
                    #Add +1 to annotation
                    count_sample[guide_no_bulge][sample][guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                    if dict_sample_to_pop[sample] not in visited_pop:
                        visited_pop.append(dict_sample_to_pop[sample])
                        count_pop[guide_no_bulge][dict_sample_to_pop[sample]][guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                    if dict_pop_to_sup[dict_sample_to_pop[sample]] not in visited_superpop:
                        visited_superpop.append(dict_pop_to_sup[dict_sample_to_pop[sample]])
                        count_superpop[guide_no_bulge][dict_pop_to_sup[dict_sample_to_pop[sample]]][guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                visited_pop = []
                visited_superpop = []
    if not found_bool:
        x[13] = 'n'
        #outFileTargets.write(line.rstrip() + '\tn\n')
    else:
        x[13] = ','.join(string_annotation)
        #outFileTargets.write(line.rstrip() + '\t' + ','.join(string_annotation) + '\n')
    last_annotation = x[13]
    last_samples = set(x[12].split(','))   
    #Save union samples + annotation
    # outFileSampleAll.write(line.rstrip() + '\t' + '\t'.join(x[-3:]) + '\n')  

    #Save cluster
    if target_scomposti_salvare:
        cluster_update.write('\t'.join(x[:-3]) + '\t' + target_scomposti_salvare[0][12] + '\t' + '\t'.join(x[13:]) + '\n')  ##This line does not contain IUPAC, needed for summary by position; Adjust sample list for target scomposed
    else:
        cluster_update.write('\t'.join(x) + '\n')       #This line does not contain IUPAC, needed for summary by position
    # cluster_update.write(line.rstrip() + '\t' + '\t'.join(x[-3:]) + '\n')   #Write line with iupac (if present)
    #Save scomposed targets
    if do_scores:
        for t in target_scomposti_salvare:
            t[13] = x[13]       #Add annotation to all scompositions
            outFileSample.write('\t'.join(t) +'\n')
            
            #Calc scores for scomposed targets
            if t[0] == 'DNA':
                cfd_score = calc_cfd(t[1][int(t[bulge_pos]):], t[2].upper()[int(t[bulge_pos]):-3], t[2].upper()[-2:], mm_scores, pam_scores)
            else:
                cfd_score = calc_cfd(t[1], t[2].upper()[:-3], t[2].upper()[-2:], mm_scores, pam_scores)
            outFileCFD.write('\t'.join(t) + '\t' + str(cfd_score) + '\n')
            sum_cfd = sum_cfd + cfd_score
            try:
                guides_dict[t[1]] = guides_dict[t[1]] + cfd_score
            except:
                guides_dict[t[1]] = cfd_score

            if t[bulge_pos + 1] == '0':    #DOENCH
                #Total = 0 -> do Classification
                for samp in t[12].split(','):
                    count_sample[guide_no_bulge][samp]['refposition'][1] += 1  #Add +1 at specific VAR onTarget
                    count_sample[guide_no_bulge][samp]['refposition'][0] += 1   #Visited +1
                #DOENCH: estraggo sequenza
                with open(outputFile + '.bedfile_tmp.bed', 'w+') as bedfile:
                    remove_tmp_bedfile = True
                    if t[6] == '+':
                        bedfile.write(t[3] + '\t' + str(int(t[4]) - 4 ) + '\t' + str(int(t[4]) + 23 + 3 ))
                    else:
                        bedfile.write(t[3] + '\t' + str(int(t[4]) - 3 ) + '\t' + str(int(t[4]) + 23 + 4 ))
                    
                extr = subprocess.Popen(['bedtools getfasta -fi ' + refgenomedir + '/' + t[3] + add_enr + add_ext + ' -bed ' + outputFile + '.bedfile_tmp.bed'], shell = True, stdout=subprocess.PIPE)  #TODO insert option for .fasta
                extr.wait()
                out, err = extr.communicate()
                out = out.decode('UTF-8')
                if t[6] == '+':
                    sequence_doench = out.strip().split('\n')[-1].upper()
                    # sequence_doench = sequence_doench[:4] + t[2] + sequence_doench[-3:]   #Uncomment to use sequence specific for sample
                else:
                    sequence_doench = reverse_complement_table(out.strip().split('\n')[-1].upper())
                    # sequence_doench = sequence_doench[:4] + t[2] + sequence_doench[-3:]   #Uncomment to use sequence specific for sample
                
                if t[1] not in targets_for_doench:
                    targets_for_doench[t[1]] = {'ref': [], 'enr': []}
                doenchForIupac(sequence_doench, t[1], 'enr')  #Get all possible targets with iupac itertools for doench
            else:
                decrease_ref_count.append(t[12]) #Save DNA, RNA, X with total != 0
            
        if not tuple_var_ref:       #Calculate scores for reference targets
            if x[0] == 'DNA':
                cfd_score = calc_cfd(x[1][int(x[bulge_pos]):], x[2].upper()[int(x[bulge_pos]):-3], x[2].upper()[-2:], mm_scores, pam_scores)
            else:
                cfd_score = calc_cfd(x[1], x[2].upper()[:-3], x[2].upper()[-2:], mm_scores, pam_scores)
            outFileCFD.write('\t'.join(x) + '\t' + str(cfd_score) + '\n')
            sum_cfd = sum_cfd + cfd_score
            try:
                guides_dict[x[1]] = guides_dict[x[1]] + cfd_score
            except:
                guides_dict[x[1]] = cfd_score

            if x[bulge_pos + 1] == '0':    #DOENCH
                #estraggo sequenza
                with open(outputFile + '.bedfile_tmp.bed', 'w+') as bedfile:
                    remove_tmp_bedfile = True
                    if x[6] == '+':
                        bedfile.write(x[3] + '\t' + str(int(x[4]) - 4 ) + '\t' + str(int(x[4]) + 23 + 3 ))
                    else:
                        bedfile.write(x[3] + '\t' + str(int(x[4]) - 3 ) + '\t' + str(int(x[4]) + 23 + 4 ))
                #Extract sequence from REFERENCE
                extr = subprocess.Popen(['bedtools getfasta -fi ' + refgenomedir + '/' + x[3] + add_enr + add_ext + ' -bed ' + outputFile + '.bedfile_tmp.bed'], shell = True, stdout=subprocess.PIPE)  #TODO insert option for .fasta
                extr.wait()
                out, err = extr.communicate()
                out = out.decode('UTF-8')
                if x[6] == '+':
                    sequence_doench = out.strip().split('\n')[-1].upper()
                    # sequence_doench = sequence_doench[:4] + x[2] + sequence_doench[-3:]
                else:
                    sequence_doench = reverse_complement_table(out.strip().split('\n')[-1].upper())
                    # sequence_doench = sequence_doench[:4] + x[2] + sequence_doench[-3:]
                
                if x[1] not in targets_for_doench:
                    targets_for_doench[x[1]] = {'ref': [], 'enr': []}
                doenchForIupac(sequence_doench, x[1], 'ref')  #Get all possible targets with iupac itertools for doench

    else:
        for t in target_scomposti_salvare:
            t[13] = x[13]       #Add annotation to all scomposition
            if t[bulge_pos + 1] == '0':
                for samp in t[12].split(','):
                    count_sample[guide_no_bulge][samp]['refposition'][1] += 1  #Add +1 at specific VAR onTarget
                    count_sample[guide_no_bulge][samp]['refposition'][0] += 1   #Visited +1
            else:
                decrease_ref_count.append(t[12])        #Save samples that have not 0 on Total column --> if in next iteration i need to save 
                                # the REF, i decrease this samples because they do not have REF of VAR On target
            outFileSample.write('\t'.join(t) + '\n')
    
    if not tuple_var_ref:
        if x[bulge_pos + 1] == '0':     #If total column is 0  -> On-target REF
            ontarget_reference_count[guide_no_bulge] +=1      #Add +1 to count ontarget for all samples (all samples have this on target)
        outFileSample.write('\t'.join(x) + '\n') #Save REF target in samples.annotation, needed for sum by guide
    decrease_ref_count = ','.join(decrease_ref_count)
    # lines_processed +=1
    # if lines_processed % (mod_tot_line) == 0:
    #     print('Annotation: Total progress ' + str(round(lines_processed /total_line *100, 2)) + '%')

############ SAVE SUMMARIES ############


#scorro tutto il dict total e scrivo il summary, targets e ogni annotation
outFileSummary.write("-Summary_Total\n")
outFileSummary.write('targets' + '\t'+'\t'.join(str(i) for i in totalDict['targets'])+'\n')
for elem in sorted(totalDict.keys(), key = lambda s : s.lower()):
    if elem == 'targets':
        continue
    outFileSummary.write(str(elem)+'\t'+'\t'.join(str(i) for i in totalDict[elem])+'\n')


for elem in guideDict.keys():
    outFileSummary.write("-Summary_"+str(elem)+'\n')
    outFileSummary.write('targets'+'\t'+'\t'.join(str(i) for i in guideDict[elem]['targets'])+'\n')
    for item in sorted(annotationsSet, key = lambda s : s.lower()):
        outFileSummary.write(str(item)+'\t'+'\t'.join(str(i) for i in guideDict[elem][item])+'\n')

#Write summaries for samples, pop, superpop
if summary_samples:
    for guide in guideDict:
        #Save sample summary
        with open(outputFile + '.sample_annotation.' + guide +'.samples.txt', 'w+') as result:
            result.write('-Summary_Total\n')
            result.write('targets'+'\t'+'\t'.join(str(i) for i in guideDict[guide]['targets'])+'\n')
            for item in sorted(annotationsSet, key = lambda s : s.lower()):
                result.write(str(item)+'\t'+'\t'.join(str(i) for i in guideDict[guide][item])+'\n')
            #Write sample specific counting, put [0]*10 if sample was not found
            for sample in all_samples:
                result.write('-Summary_' + sample + '\n')
                try:
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in count_sample[guide][sample]['targets']) + '\n')
                except: #Sample not found in targets
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')
                for item in sorted(annotationsSet, key = lambda s : s.lower()):
                    try:
                        result.write(item + '\t' + '\t'.join(str(i) for i in count_sample[guide][sample][item]) + '\n')
                    except:
                        result.write(item + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')
        
        #Save population summary
        with open(outputFile + '.sample_annotation.' + guide +'.population.txt', 'w+') as result:
            result.write('-Summary_Total\n')
            result.write('targets'+'\t'+'\t'.join(str(i) for i in guideDict[guide]['targets'])+'\n')
            for item in sorted(annotationsSet, key = lambda s : s.lower()):
                result.write(str(item)+'\t'+'\t'.join(str(i) for i in guideDict[guide][item])+'\n')
            #Write population specific counting, put [0]*10 if sample was not found
            for population in set(all_pop):
                result.write('-Summary_' + population + '\n')
                try:
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in count_pop[guide][population]['targets']) + '\n')
                except: #Sample not found in targets
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')
                for item in sorted(annotationsSet, key = lambda s : s.lower()):
                    try:
                        result.write(item + '\t' + '\t'.join(str(i) for i in count_pop[guide][population][item]) + '\n')
                    except:
                        result.write(item + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')
        
        #Save superpopulation summary
        with open(outputFile + '.sample_annotation.' + guide +'.superpopulation.txt', 'w+') as result:
            result.write('-Summary_Total\n')
            result.write('targets'+'\t'+'\t'.join(str(i) for i in guideDict[guide]['targets'])+'\n')
            for item in sorted(annotationsSet, key = lambda s : s.lower()):
                result.write(str(item)+'\t'+'\t'.join(str(i) for i in guideDict[guide][item])+'\n')
            #Write superpopulation specific counting, put [0]*10 if sample was not found
            for superpop in superpopulation:
                result.write('-Summary_' + superpop + '\n')
                try:
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in count_superpop[guide][superpop]['targets']) + '\n')
                except: #Sample not found in targets
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')
                for item in sorted(annotationsSet, key = lambda s : s.lower()):
                    try:
                        result.write(item + '\t' + '\t'.join(str(i) for i in count_superpop[guide][superpop][item]) + '\n')
                    except:
                        result.write(item + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')
    with open(outputFile + '.PopulationDistribution.txt', 'w+') as pop_distribution:
        for g in guideDict:
            pop_distribution.write('-Summary_' + g + '\n')
            pop_distribution.write('REFERENCE' + '\t' + '\t'.join([ ','.join(str(v) for v in t) for t in add_to_general_table[g]['distributions']]) + '\n')
            for superpop in superpopulation:
                try:
                    pop_distribution.write(superpop + '\t' + '\t'.join([ ','.join(str(v) for v in t) for t in count_superpop[g][superpop]['distributions']]) + '\n')
                except:
                    pop_distribution.write(superpop + '\t' + '\t'.join([ ','.join(str(v) for v in t) for t in [[0] * (max_bulges + 1) for x in range(10)]]) + '\n')

                #[[0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0]] ->   0,0,0   0,0,0   0,0,0   0,0,0   0,0,0 
#Write sumref for barplot for targets in top1 form of var/ref search
if summary_barplot_from_total:
    with open(outputFile + '.sumref.Annotation.summary.txt', 'w+') as result:
        result.write('-Summary_Total\n')
        result.write('targets'+'\t'+'\t'.join(str(i - count_unique['targets'][pos]) for pos,i in enumerate(totalDict['targets'])) + '\n')
        for elem in sorted(annotationsSet, key = lambda s : s.lower()):
            result.write(str(elem)+'\t'+'\t'.join(str(i - count_unique[elem][pos]) for pos, i in enumerate(totalDict[elem]))+'\n')
        for guide in count_unique_for_guide:
            result.write('-Summary_' + guide + '\n')
            result.write('targets' + '\t' + '\t'.join(str(i - count_unique_for_guide[guide]['targets'][pos]) for pos,i in enumerate(guideDict[guide]['targets'])) + '\n')
            for annotation in sorted(annotationsSet, key = lambda s : s.lower()):
                result.write(annotation + '\t' + '\t'.join(str(i - count_unique_for_guide[guide][annotation][pos]) for pos, i in enumerate(guideDict[guide][annotation])) + '\n')

#SAVE SCORES#
all_scores = []
if do_scores:
    with open( outputFile + '.scores.txt', 'w+') as res, open(sys.argv[8], 'r') as guides:
        man = multiprocessing.Manager()
        shared_doench = man.list() #list containing max doech for each thread
        guides = guides.read().strip().split('\n')
        for g in guides:
            guides_dict_doench[g] = {'ref': 0, 'enr': 0} 
            if g not in guides_dict:
                guides_dict[g] = 0    
            if g not in targets_for_doench:
                guides_dict_doench[g] = {'ref': 0, 'enr': 0} 
            else:
                for gentype in ['ref', 'enr']:
                    if len(targets_for_doench[g][gentype]) == 0:
                        guides_dict_doench[g][gentype] = 0
                        continue
                    if len (targets_for_doench[g][gentype]) > SIZE_DOENCH:
                        jobs = []
                        remaining_splits = (len(targets_for_doench[g][gentype])//SIZE_DOENCH) + 1
                        for i in range ((len(targets_for_doench[g][gentype])//SIZE_DOENCH) + 1):
                            for thr in range (min(N_THR, remaining_splits)):
                                p = multiprocessing.Process(target = doenchParallel, args=(np.asarray(targets_for_doench[g][gentype][i*N_THR*SIZE_DOENCH + thr*SIZE_DOENCH : min( i*N_THR*SIZE_DOENCH + (thr+1)*SIZE_DOENCH,len(targets_for_doench[g][gentype]))]), model, shared_doench,) )
                                remaining_splits -= 1
                                p.start()
                                jobs.append(p)
                            for i in jobs:
                                i.join()
                        
                        guides_dict_doench[g][gentype] = max(shared_doench)
                        shared_doench =  man.list()
                    else:
                        start_time = time.time()
                        doench_score =  azimuth.model_comparison.predict(np.asarray(targets_for_doench[g][gentype]), None, None, model= model, pam_audit=False)
                        doench_score = [np.around(i * 100) for i in doench_score]
                        guides_dict_doench[g][gentype] =  int(max(doench_score))
            if guides_dict[g] == 0:     #NO CFD CALCULATED
                all_scores.append([g, 0, str(guides_dict_doench[g]['ref']), str(max(guides_dict_doench[g]['ref'], guides_dict_doench[g]['enr']))])
            else:
                all_scores.append([g, int(round((100/(100 + guides_dict[g]))*100)), str(guides_dict_doench[g]['ref']), str(max(guides_dict_doench[g]['ref'], guides_dict_doench[g]['enr'])) ])
        all_scores.sort(key = lambda x: x[1], reverse = True)
        res.write('#crRNA\tCFD\tDoench 2016 - Reference\tDoench 2016 - Enriched\n')
        for guide_scored in all_scores:  
            res.write('\t'.join([str(el) for el in guide_scored]) + '\n') 

#Save additional values from semicommon for general guide table
with open(outputFile + '.addToGeneralTable.txt', 'w+') as add_file:
    for guide in add_to_general_table:
        add_file.write(guide + '\t' + '\t'.join([str(x) for x in add_to_general_table[guide]['add']]) + '\n')

#Calculate and save classes for each sample
total_class = dict()    #For each guide, Save total of class [0, 1 , 1+]
with open(outputFile + '.SampleClasses.txt', 'w+') as sample_class:
    #Header
    sample_class.write('Sample')
    for g in count_sample.keys():
        total_class[g] = [0, 0, 0, 0]   #0 0+ 1 1+
        sample_class.write('\t' + g)
    sample_class.write('\n')
    #Write classes
    for s in all_samples:
        sample_class.write(s)
        for g in count_sample.keys():
            # count_on_target = (count_sample[g][s]['refposition'][1] + count_sample[g][s]['refposition'][2] + ontarget_reference_count[g]) #Sum created new ontarget + destroyed ontarget + ontarget on all samples
            if s not in count_sample[g]:
                if ontarget_reference_count[g] > 0:
                    sample_class.write('\t1')
                    total_class[g][2] += 1
                else:
                    sample_class.write('\t0')   #If there are no REF, the samples has 0 on targets
                    total_class[g][0] += 1
                continue
            if count_sample[g][s]['refposition'][1] > 0: #CLASS 1+ : sample have at least 1 new on target w.r.t the reference 
                sample_class.write('\t1+')
                total_class[g][3] += 1
            elif count_sample[g][s]['refposition'][0] == 0 : #CLASS 1: sample have the same ontarget as ref
                if ontarget_reference_count[g] > 0:
                    sample_class.write('\t1')
                    total_class[g][2] += 1
                else:
                    sample_class.write('\t0')   #If there are no REF, the samples has 0 on targets
                    total_class[g][0] += 1
            elif ontarget_reference_count[g] + count_sample[g][s]['refposition'][2] == 0: #Class 0: all the reference targets were destroyed
                sample_class.write('\t0')       #NOTE samples that have at least one new OnTarget, but ontarget_ref + ['refpos'][2] == 0 are counted as 1+
                total_class[g][0] +=1
            else:
                sample_class.write('\t0+')      #Class 0+: only some REF were destroyed, sample is still targetable
                total_class[g][1] +=1
        sample_class.write('\n')
    sample_class.write('Total for Class 0 - 0+ - 1 - 1+')
    for g in count_sample.keys():
        sample_class.write('\t' + '-'.join([str(x) for x in total_class[g]]))
if total_error > 0:
    print('Skipped SNP:', total_error)
if remove_tmp_bedfile:
    os.remove(outputFile + '.bedfile_tmp.bed')
if 'Step' not in step:
    print("SAMPLE PROCESSING AND ANNOTATION COMPLETED IN: %s seconds" % (time.time() - start_time_total))

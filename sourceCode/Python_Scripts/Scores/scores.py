#!/usr/bin/env python

# Faster script + multiprocessing
# Script that calculates cfd score for targets with bulge 'X', and save the accumulated cfd score for the input guide
# Also calculates the Doench score if the targets has bulge 'X' and 0 mms (doench = 0 if no such target exists)

#argv 1 = target file
#argv 2 is genome_directory (eg ../../Genomes/hg19/)
#argv 3 is pam file -> to check if len is 23 and pam is NGG
#argv 4 is guide file
import time
import pickle
import re
import sys
import os
from os import listdir
from os.path import isfile, join
import numpy as np
import subprocess
import azimuth.model_comparison
import string
import itertools
import multiprocessing
SIZE_DOENCH = 10000
N_THR = 3
# doench_string.append(seq)
# doench_score =  azimuth.model_comparison.predict(np.asarray(doench_string), None, None, model= model, pam_audit=False)
# doench_score = np.around(doench_score * 100)

def doenchParallel(targets, model, result):
  start_time = time.time()
  doench_score =  azimuth.model_comparison.predict(targets,None, None, model= model, pam_audit=False)
  doench_score = [np.around(i * 100) for i in doench_score]
  max_doench = int(max(doench_score))
  result.append(max_doench)


def doenchForIupac(sequence_doench, guide_seq):
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
        targets_for_doench[guide_seq].append(''.join(t))
  else:
    targets_for_doench[guide_seq].append(sequence_doench)


def get_mm_pam_scores():
  try:
    mm_scores = pickle.load(open(os.path.dirname(os.path.realpath(__file__)) + '/mismatch_score.pkl', 'rb'))
    pam_scores = pickle.load(open(os.path.dirname(os.path.realpath(__file__)) +'/PAM_scores.pkl', 'rb'))
    return (mm_scores, pam_scores)
  except:
    raise Exception("Could not find file with mismatch scores or PAM scores")


def revcom(s):
  basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
  letters = list(s[::-1])
  letters = [basecomp[base] for base in letters]
  return ''.join(letters)


# Calculates CFD score
def calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores):
    
    score = 1
    dna_gp = 0
    sg = sg.replace('T', 'U')
    guide_seq = guide_seq.replace('T', 'U')
    s_list = list(sg)
    guide_seq_list = list(guide_seq)
    for i, sl in enumerate(s_list):
      
      if guide_seq_list[i] == sl:

          score *= 1

      else:
          key = 'r' + guide_seq_list[i] + ':d' + revcom(sl) + ',' + str(i + 1)
          score *= mm_scores[key]
          if '-' in guide_seq_list[i]:
            dna_gp = dna_gp + 1
      
    score *= pam_scores[pam]
    
    return score

tab = str.maketrans("ACTGRYSWMKHDBVactgryswmkhdbv", "TGACYRSWKMDHVBtgacyrswkmdhvb") 

def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]
#if __name__ == '__main__':

mm_scores, pam_scores = get_mm_pam_scores()
guides_dict = dict()
guides_dict_doench = dict()
targets_for_doench = dict()

score_filename = sys.argv[1].strip().split('.targets.txt')[0] + '.scores.txt'
bedfile_tmp_name = sys.argv[1].strip().split('.targets.txt')[0] + '.bedfile_tmp.bed'
remove_bedfile_tmp = False
with open(sys.argv[3]) as pamfile:
  line = pamfile.readline().strip().split(' ')
  if len(line[0]) != 23 or 'NGG' not in line[0]:
    with open(score_filename, 'w+') as result:
      result.write('NO SCORES')
      exit()

N_THR = multiprocessing.cpu_count() // 2

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

start = time.time()

#Get list of chromosomes file to set enr_str
chromosome_files = [f for f in listdir(sys.argv[2]) if isfile(join(sys.argv[2], f))]
# enr = sys.argv[2].split('/')
enr_str = ''
chr_ext = '.fa'
if '.enriched.' in chromosome_files[0]:
  enr_str = '.enriched'
elif '.indels.' in chromosome_files[0]:
  enr_str = '.indels'
if '.fasta' in chromosome_files[0]:
  chr_ext = '.fasta'
# if enr[-1]:
#   if'+' in enr[-1]:
#     enr_str = '.enriched'
# else:
#   if'+' in enr[-2]:
#     enr_str = '.enriched'

with open( os.path.dirname(os.path.realpath(__file__)) + "/azimuth/saved_models/V3_model_nopos.pickle", 'rb') as f:
  model = pickle.load(f)
max_doench = 0
n_of_acceptable_cfd = 0
sum_cfd = 0
cfd_scores = []

# NOTE uncomment for progress bas
# process = subprocess.Popen(['wc', '-l', sys.argv[1]], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# out, err = process.communicate()
# total_line = int(out.decode('UTF-8').split(' ')[0])
# if total_line < 2:
#     print('WARNING! Input file has no targets')
#     sys.exit()
# if total_line < 10:
#     mod_tot_line = 1
# else:
#     mod_tot_line = int(total_line/90)
# lines_processed = 0


all_word = []
with open (sys.argv[1]) as result:
  
  #Calc CDF score
  for target  in result:
    target = target.strip().split('\t')
    if 'X' not in target[0]:
      continue
    
    guide_seq = target[1]
    off = target[2].upper()
    m_guide_seq = re.search('[^ATCGN-]', guide_seq)
    m_off = re.search('[^ATCG-]', off)  
    iup_off = []
    first = True
    start_iup_off = 1
    
    if (m_guide_seq is None) and (m_off is None):
       
      #Calc CFD
          
      
      pam = off[-2:]  
      sg = off[:-3]
      #print("off. ", off)
      #print ("sg: ", sg)
      #print ("guide_seq: ", guide_seq)
      
      cfd_score = calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores)
      if (target[7] == '0'):    #TODO se cambio inserendo pos cluister, devo cambiareanche qui, da 6 a 7 (con colonna pos cluster)
        #estraggo sequenza
        with open(bedfile_tmp_name, 'w+') as bedfile:
          remove_bedfile_tmp = True
          if target[6] == '+':
            bedfile.write(target[3] + '\t' + str(int(target[4]) - 4 ) + '\t' + str(int(target[4]) + 23 + 3 ))
          else:
            bedfile.write(target[3] + '\t' + str(int(target[4]) - 3 ) + '\t' + str(int(target[4]) + 23 + 4 ))
        
        extr = subprocess.Popen(['bedtools getfasta -fi ' + sys.argv[2] + '/' + target[3] +  enr_str + chr_ext + ' -bed ' + bedfile_tmp_name], shell = True, stdout=subprocess.PIPE) 
        extr.wait()
        out, err = extr.communicate()
        out = out.decode('UTF-8')
        if target[6] == '+':
          sequence_doench = out.strip().split('\n')[-1].upper()
        else:
          sequence_doench = reverse_complement_table(out.strip().split('\n')[-1].upper())
        
        if target[1] not in targets_for_doench:
          targets_for_doench[target[1]] = []
        doenchForIupac(sequence_doench, target[1])  #Get all possible targets with iupac itertools for doench

      sum_cfd = sum_cfd + cfd_score
      try:
        guides_dict[target[1]] = guides_dict[target[1]] + cfd_score
      except:
        guides_dict[target[1]] = cfd_score
      if cfd_score > 0.023:
        n_of_acceptable_cfd = n_of_acceptable_cfd +1  
      
    else:
      if "N" in off:
        continue
      if (target[7] == '0'):  #NOTE change from 6 to 7 if input file has cluster position column
        with open(bedfile_tmp_name, 'w+') as bedfile:
          remove_bedfile_tmp = True
          if target[6] == '+':
            bedfile.write(target[3] + '\t' + str(int(target[4]) - 4 ) + '\t' + str(int(target[4]) + 23 + 3 ))
          else:
            bedfile.write(target[3] + '\t' + str(int(target[4]) - 3 ) + '\t' + str(int(target[4]) + 23 + 4 ))
        
        extr = subprocess.Popen(['bedtools getfasta -fi ' + sys.argv[2] + '/' + target[3] + enr_str + chr_ext + ' -bed ' + bedfile_tmp_name], shell = True, stdout=subprocess.PIPE) 
        extr.wait()
        out, err = extr.communicate()
        out = out.decode('UTF-8')
        if target[6] == '+':
          sequence_doench = out.strip().split('\n')[-1].upper()
        else:
          sequence_doench = reverse_complement_table(out.strip().split('\n')[-1].upper())
        
        if target[1] not in targets_for_doench:
          targets_for_doench[target[1]] = []
        doenchForIupac(sequence_doench, target[1])  #Get all possible targets with iupac itertools for doench

      i = 0
      for char in off:
        if char in iupac_code:
          n = len(iup_off)
          for list_char in iupac_code[char]:
            if not first:  
              for test in range(n - start_iup_off, n):
                iup_off.append(iup_off[test][:i] + list_char + iup_off[test][i+1:])
                
              
            else:
              iup_off.append(off[:i] + list_char + off[i+1:])
          first = False
          start_iup_off = start_iup_off * len(iupac_code[char])
          
        i += 1
      dna_gap_removal = True
      for no_iup_str in range(len(iup_off) - start_iup_off, len(iup_off)):
        
        
        no_iup_gap_srt = iup_off[no_iup_str] #se non ci sono gap passo la stringa non modificata al calc_cfd
      
        #Calc CFD
    

        pam = no_iup_gap_srt[-2:]   
        sg = no_iup_gap_srt[:-3]
        
        cfd_score = calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores)
        sum_cfd = sum_cfd + cfd_score
        try:
          guides_dict[target[1]] = guides_dict[target[1]] + cfd_score
        except:
          guides_dict[target[1]] = cfd_score
    
    #NOTE decomment for progress bar
    # lines_processed +=1
    # if lines_processed % (mod_tot_line) == 0:
    #     print('Scoring: Total progress ' + str(round(lines_processed /total_line *100, 2)) + '%')

job_id = sys.argv[1].split('/')[-1].split('.')[0]

all_scores = []
with open( score_filename, 'w+') as res, open(sys.argv[4], 'r') as guides:
  man = multiprocessing.Manager()
  shared_doench = man.list() #list containing max doech for each thread
  guides = guides.read().strip().split('\n')
  for g in guides:
    guides_dict_doench[g] = 0
    if g not in guides_dict:
      guides_dict[g] = 0    
  #for k in guides_dict.keys():
    if g not in targets_for_doench:
      guides_dict_doench[g] = 0
    else:
      if len (targets_for_doench[g]) > SIZE_DOENCH:
        jobs = []
        remaining_splits = (len(targets_for_doench[g])//SIZE_DOENCH) + 1
        for i in range ((len(targets_for_doench[g])//SIZE_DOENCH) + 1):
          for thr in range (min(N_THR, remaining_splits)):
            p = multiprocessing.Process(target = doenchParallel, args=(np.asarray(targets_for_doench[g][i*N_THR*SIZE_DOENCH + thr*SIZE_DOENCH : min( i*N_THR*SIZE_DOENCH + (thr+1)*SIZE_DOENCH,len(targets_for_doench[g]))]), model, shared_doench,) )
            remaining_splits -= 1
            p.start()
            jobs.append(p)
          for i in jobs:
            i.join()
          
        guides_dict_doench[g] = max(shared_doench)
        shared_doench =  man.list()
      else:
        start_time = time.time()
        doench_score =  azimuth.model_comparison.predict(np.asarray(targets_for_doench[g]), None, None, model= model, pam_audit=False)
        doench_score = [np.around(i * 100) for i in doench_score]
        guides_dict_doench[g] =  int(max(doench_score))
    if guides_dict[g] == 0:   #NO CFD WAS CALCULATED
      all_scores.append([g, 0, str(guides_dict_doench[g]) ])
    else:  
      all_scores.append([g, int(round((100/(100 + guides_dict[g]))*100)), str(guides_dict_doench[g]) ])
  all_scores.sort(key = lambda x: x[1], reverse = True)
  res.write('#crRNA\tCFD\tDoench 2016\n')
  for guide_scored in all_scores:  
    res.write('\t'.join([str(el) for el in guide_scored]) + '\n') 

if remove_bedfile_tmp:
  os.remove(bedfile_tmp_name)
end = time.time()

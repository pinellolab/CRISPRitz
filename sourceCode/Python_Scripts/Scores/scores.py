#!/usr/bin/env python

import time
import pickle
import re
import sys
import os
import azimuth.model_comparison
import numpy as np
import subprocess

DOENCH_SIZE = 10000
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


#if __name__ == '__main__':

#argv 1 = directory containing fasta files (input by user), absolute path
#argv 2 = pam+guide length
#argv 3 = target filename, absolute path
#argv 4 = true if pam at start. 

fasta_files = [f for f in os.listdir(sys.argv[1]) if os.path.isfile(os.path.join(sys.argv[1], f))]

origin_path = os.path.dirname(os.path.realpath(__file__))

target_name = sys.argv[3].split("/")[-1].split(".targets")[0]
if (sys.argv[4] == "True"):
  pam_at_start = True
else:
  pam_at_start = False

#Check if unwanted files (.bed, .ordered, .seqs) are in the directory
if (os.path.isfile(origin_path + "/" + target_name + ".order.txt")):
  os.remove(origin_path + "/" + target_name + ".order.txt")
if (os.path.isfile(origin_path + "/" + target_name + ".bed")):
  os.remove(origin_path + "/" + target_name + ".bed")
if (os.path.isfile(origin_path + "/" + target_name + ".seqs.txt")):
  os.remove(origin_path + "/" + target_name + ".seqs.txt")


#Check if target file has at least 1 target
if (int(subprocess.run(["wc", "-l" , sys.argv[3]], stdout=subprocess.PIPE).stdout.decode('utf-8').split()[0]) == 1):
  print("No targets were found. Skipping score calculation")
  sys.exit()  
 


print("Preprocessing data...")
#Sort target file and create an ordered by chr .bed 
start_sort = time.time()
subprocess.run(["sort", "-k4,4", sys.argv[3], "-o" , origin_path + "/" + target_name + ".order.txt"]) #salvato in CDF
end_sort = time.time()
#print("Sort: ", end_sort - start_sort)

subprocess.run(["sed -i '/^#/ d' " + origin_path + "/" + target_name + ".order.txt"], shell=True)


start_getseq = time.time()
subprocess.run([origin_path+"/get_sequences.sh", origin_path + "/" + target_name + ".order.txt", sys.argv[2], origin_path + "/" + target_name, sys.argv[4]])  #.bed saved in Scores directory          
end_getseq = time.time()
#print("Get seq: ", end_getseq - start_getseq)

#Order fasta_file and do bamtools getfasta
start_bedtools = time.time()
fasta_files.sort()
for fasta in fasta_files:
  if fasta.endswith(".fasta") or fasta.endswith(".fa") :
    subprocess.run(["bedtools getfasta -fi " + sys.argv[1] + fasta + " -bed " + origin_path + "/" + target_name + ".bed" + " >> " + origin_path + "/" + target_name + ".seqs.txt"], shell= True, stderr=subprocess.DEVNULL)

if (int(subprocess.run(["wc", "-l" , origin_path + "/" + target_name + ".seqs.txt"], stdout=subprocess.PIPE).stdout.decode('utf-8').split()[0]) == 0):
  print("WARNING! No sequences extracted, check if chromosomes names are the same in the fasta file and in the targets file")
  sys.exit()
end_bedtools = time.time()
#print("Bedtools: ", end_bedtools - start_bedtools)

#Remove indices created by bedtools
fasta_files = [f for f in os.listdir(sys.argv[1]) if os.path.isfile(os.path.join(sys.argv[1], f))]
for fai in fasta_files:
  if fai.endswith(".fai"):
    os.remove(sys.argv[1] + fai)    

start_sed = time.time()
subprocess.run(["sed -i '1~2d' " + origin_path + "/" + target_name + ".seqs.txt"], shell = True) 
end_sed = time.time()

#print("Sed time: ", end_sed - start_sed)

#Get number of lines of target file, to have a % of scores calculated
num_lines = int(subprocess.run(["wc", "-l" , sys.argv[3]], stdout=subprocess.PIPE).stdout.decode('utf-8').split()[0])
steps_end_calc = 0
total_steps = np.ceil(num_lines / DOENCH_SIZE) #DOENCH_SIZE is size of batch for Doench score

result = open( sys.argv[3].split('.targets')[0] + ".scores.txt",mode= "w")  #create a file with CDF score with name based on input file
target_len = int(sys.argv[2])
mm_scores, pam_scores = get_mm_pam_scores()


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
          "V":("A", "C", "G")
          }

start = time.time()

result.write("#Bulge type\tcrRNA\tDNA\tChromosome	Position\tDirection\tMismatches\tBulge Size\tCFD\tDoench 2016\n")
n_of_acceptable_cfd = 0

cfd_scores = []
doench_string = []
# doench_string_dna_1 = []
# doench_string_dna_2 = []
# doench_string_rna_1 = []
# doench_string_rna_2 = []
size_doe = 0
all_word = []
with open (origin_path + "/" + target_name + ".seqs.txt") as file_seqs, open(origin_path + "/" + target_name + '.order.txt') as file_order, open(origin_path + "/azimuth/saved_models/V3_model_nopos.pickle", 'rb') as f:
  model = pickle.load(f)
  #Calc CDF score
  for seq  in file_seqs:
    word=file_order.readline().split("\t")
    word[-1] = word[-1].strip()
    seq = seq.upper().strip()
    guide_seq = word[1]
    off = seq.upper()
    m_guide_seq = re.search('[^ATCGN-]', guide_seq)
    m_off = re.search('[^ATCG-]', off)  
    iup_off = []
    first = True
    start_iup_off = 1
    
    if (m_guide_seq is None) and (m_off is None):
      size_doe = size_doe + 1
      # gap_pos = [None]
      # if "DNA" in word[0]:    #gap in guide
      #   gap_pos = [i for i, letter in enumerate(word[1]) if letter == "-"]
      #   guide_seq = guide_seq.replace('-', '')
      #   for g in gap_pos:
      #     off = off[0:g] + off[g+1:]
      # if "RNA" in word[0]:    #gap in target
      #   gap_pos = [i for i, letter in enumerate(word[2]) if letter == "-"]
      #   #off = off.replace('-', '')
      #   for g in gap_pos:
      #     guide_seq = guide_seq[0:g] + guide_seq[g+1:]
      #     off = off[0:g] + off[g + 1 :]
      
      #Calc CFD
      if pam_at_start:
        off = off[::-1]
        seq = seq[::-1]
      
      
      pam = off[-2:]  
      sg = off[4:-6]
      #print("off. ", off)
      #print ("sg: ", sg)
      #print ("guide_seq: ", guide_seq)
      
      cfd_score = calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores)
      # print (cfd_score)
      
      if cfd_score > 0.023:
        n_of_acceptable_cfd = n_of_acceptable_cfd +1
      
      all_word.append(word)
      cfd_scores.append(cfd_score)
      doench_string.append(seq)    
      
    else:
      if "N" in off:
        continue
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
        
        size_doe = size_doe + 1
        # gap_pos = [None]
        no_iup_gap_srt = iup_off[no_iup_str] #se non ci sono gap passo la stringa non modificata al calc_cfd
        # if "DNA" in word[0]:    #gap in guide
        #   gap_pos = [i for i, letter in enumerate(word[1]) if letter == "-"]
        #   guide_seq = guide_seq.replace('-', '')
        #   for g in gap_pos:
        #     no_iup_gap_srt = iup_off[no_iup_str][0:g] + iup_off[no_iup_str][g+1:]

        # if "RNA" in word[0]:    #gap in target
        #   gap_pos = [i for i, letter in enumerate(word[2]) if letter == "-"]
        #   for g in gap_pos:
        #     if dna_gap_removal:
        #       guide_seq = guide_seq[0:g] + guide_seq[g+1:]
        #       dna_gap_removal = False
        #     no_iup_gap_srt = iup_off[no_iup_str][0:g] + iup_off[no_iup_str][g+1:]
        
        #Calc CFD
        if pam_at_start:
          off = no_iup_gap_srt[::-1]
          iup_off[no_iup_str] = iup_off[no_iup_str][::-1]

        pam = no_iup_gap_srt[-2:]   
        sg = no_iup_gap_srt[4:-6]
        
        cfd_score = calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores)

        temp_word = word.copy()
        cfd_scores.append(cfd_score)
        all_word.append(temp_word)
        doench_string.append(iup_off[no_iup_str])


        #Restore lowercase letters
        for count, c in enumerate(all_word[-1][2], 4):
          if c.islower():
            iup_off[no_iup_str] = iup_off[no_iup_str][:count] + iup_off[no_iup_str][count].lower() + iup_off[no_iup_str][count + 1 :]  
          

        all_word[-1][2] = iup_off[no_iup_str][4:-3]
        

    if size_doe >= DOENCH_SIZE:
      start_doench = time.time()
      doench_score =  azimuth.model_comparison.predict(np.asarray(doench_string), None, None, model= model, pam_audit=False)
      doench_score = np.around(doench_score * 100)
      end_doench = time.time()
      #print("Doench: ", end_doench - start_doench)
      steps_end_calc = steps_end_calc + 1
      print("Calculating scores: " + str(np.around((steps_end_calc/total_steps) * 100, decimals= 1)) + "%")

      #save to file
      for kk in range (size_doe):
        result.write(all_word[kk][0] + "\t" + all_word[kk][1] + "\t" + all_word[kk][2] + "\t" + all_word[kk][3 ] + "\t" + all_word[kk][4] + "\t" 
          + all_word[kk][5] + "\t" + all_word[kk][6] + "\t" + all_word[kk][7] + "\t" + str(cfd_scores[kk]) + "\t" + str(doench_score[kk]) + "\n" )
      size_doe = 0
      cfd_scores = []
      all_word = []
      doench_string = []
  if doench_string:
  #eof, calculate doench for remaining string in all the lists
    doench_score =  azimuth.model_comparison.predict(np.asarray(doench_string), None, None, model= model, pam_audit=False)
    doench_score = np.around(doench_score * 100)
    print("Calculating scores: " +  "100%")

    #save to file
    for kk in range (size_doe):
      result.write(all_word[kk][0] + "\t" + all_word[kk][1] + "\t" + all_word[kk][2] + "\t" + all_word[kk][3 ] + "\t" + all_word[kk][4] + "\t" 
        + all_word[kk][5] + "\t" + all_word[kk][6] + "\t" + all_word[kk][7] + "\t" + str(cfd_scores[kk]) + "\t" + str(doench_score[kk]) + "\n" )
    cfd_scores = []
    all_word = []
    doench_string = []
  else:
    print("Calculating scores: " +  "100%")



end = time.time()
#print ("CFD score time: ", end-start)
#print ("CFD acceptable: ", n_of_acceptable_cfd)

result.close()
#Remove tmp file
os.remove(origin_path + "/" + target_name + ".order.txt")
os.remove(origin_path + "/" + target_name + ".bed")
os.remove(origin_path + "/" + target_name + ".seqs.txt")

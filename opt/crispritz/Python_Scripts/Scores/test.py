import azimuth.model_comparison
import numpy as np
import pickle
import random 
import string
import time
import os

## Generate a random string of specific characters 
def randString(length=5):
    #put your letters in the following string
    your_letters='ATCG'
    return ''.join((random.choice(your_letters) for i in range(length)))

origin_path = os.path.dirname(os.path.realpath(__file__))

sequences = np.array(['AAAAAAAAAAAAAAAAAAATGCAGCGGGAG', 'AAAAAAAAAAAAAAGACTCACAGGTGGACC', 'AAAAAAAAAAAAAATAGATATAATTGGCCT', 'CAGCTGATCTCCAGATATGACCATGGGTTT', 'CCAGAAGTTTGAGCCACAAACCCATGGTCA'
            , 'TTTGTCCCCAATCTGGGCGCGCGCCGGCGC', 'TGGCGGCCTAAGGACTCGGCGCGCCGGAAG', 'CCTAAGGACTCGGCGCGCCGGAAGTGGCCA', 'GGCGACCTCGGCTCACAGCGCGCCCGGCTA', 'TATCGCCGCGCTCGTCGTCGACAACGGCTC',
            'GGACTCGGCGCGCCGGAAGTGGCCAGGGCG', 'GACTCGGCGCGCCGGAAGTGGCCAGGGCGG', 'GCGCCCCCTGGCGGCCTAAGGACTCGGCGC', 'CCCGGCTATTCTCGCAGCTCACCATGGATG', 'GCGCCGGCGCCCCCTGGCGGCCTAAGGACT',
             'TCGGCGCGCCGGAAGTGGCCAGGGCGGGGG', 'GGCGCGCCGGAAGTGGCCAGGGCGGGGGCG', 'CCGGCATGTGCAAGGCCGGCTTCGCGGGCG', 'CGGGCCGTCTTCCCCTCCATCGTGGGGCGC',
             'GGGGCGCCGGCGCGCGCCCAGATTGGGGAC', 'CCACTTCCGGCGCGCCGAGTCCTTAGGCCG', 'GCGCCGAGTCCTTAGGCCGCCAGGGGGCGC', 'CGGAGCCGTTGTCGACGACGAGCGCGGCGA',
             'GGGGGCGCCGGCGCGCGCCCAGATTGGGGA', 'AGTCCTTAGGCCGCCAGGGGGCGCCGGCGC', 'ACGATGGAGGGGAAGACGGCCCGGGGGGCA',
             'TTCTACTACTACAGTTGAGTGGTCTGGCAT', 'AAAAAAGAGTTCCAGAGTGCTCTGGGGTTT', 'AAAAAAAAGAGTTCCAGAGTGCTCTGGGGT', 'TTCTCTAGAATGCAACATTTTTGTAGGACA'])
# sequences = np.array([randString(29))
# for i in range(100):
#     ins = randString(29)
#     sequences = np.append(sequences,ins)     
#     if (i%1000==0):
#         print (i)
print("DONE")        

#amino_acid_cut_positions = np.array([2, 2, 4])
#percent_peptides = np.array([0.18, 0.18, 0.35])
amino_acid_cut_positions = None
percent_peptides = None

res = np.array([0,0,0,0,0,0,0,47,43,57,40,52,33,55,29,43, 23, 61,64, 22, 47, 65, 62, 20, 34, 52, 49,75, 50, 27])

with open(origin_path + "/azimuth/saved_models/V3_model_nopos.pickle", 'rb') as f:
            model = pickle.load(f)
start = time.time()
predictions = azimuth.model_comparison.predict(sequences, amino_acid_cut_positions, percent_peptides, model=model, pam_audit=False)
end = time.time()
print (end-start)
for i, prediction in enumerate(predictions):
    print (sequences[i], prediction, res[i])
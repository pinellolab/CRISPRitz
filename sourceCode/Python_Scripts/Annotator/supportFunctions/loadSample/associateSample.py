'''
Script that loads the association
SAMPLE -> POPULATION -> SUPERPOPULATION
and
SUPERPOPULATION -> LIST OF POPULATION ; POPULATION -> LIST OF SAMPLE
Also return list of all samples, all population, all superpopulations
for the given input file. The input file is a txt file with four columns, separated by a tabulation, where the first column
represents the sample ID, the second column the associated Population ID,  the third column the SUuperpopulation ID, and the last column the 
sample Gender (if last column is not available, 'n/a' is written). 
The script returns a dictionary sample_to_pop, pop_to_superpop, superpop_to_pop, pop_to_sample, and the lists
all_samples, all_pop, all_superpop and list_gender or 
None * 8 if input file does not exists or is in incorrect format
'''
import sys
import os
from os.path import isfile, isdir,join      

#argv 1 is input .txt, consisting of 3 columns for sample, population and superpopulation

def loadSampleAssociation(id_sample_file):
    '''
    Given in input a file, it checks for structure correctness and returns 4 dictionaries, 4 lists, 1 dictionary:
    - SAMPLE -> POPULATION
    - POPULATION -> SUPERPOPULATION
    - SUPERPOPULATION -> LIST OF POPULATIONS
    - POPULATIONS -> LIST OF SAMPLES
    - ALL SAMPLES
    - ALL POPULATIONS
    - ALL SUPERPOPULATIONS
    - SAMPLE -> GENDER
    '''
    if not isfile(id_sample_file):
        print('Warning! The sample ID file does not exists')
        print('Exit..')
        sys.exit() 
        return None, None

    sample_to_pop = dict()
    pop_to_superpop = dict() 
    superpop_to_pop = dict()   #For each superpopulation returns list of population
    pop_to_sample = dict()      #For each population return list of sample
    all_samples = set()
    all_pop= set()
    all_superpop = set()
    gender_sample = dict()
    with open (id_sample_file) as in_file:
        #Check correct format of file
        line = next(in_file)
        if '#' in line:
            line = next(in_file)    #Skip header
        
        line = line.strip().split('\t')

        if len(line) < 3:
            print('Warning! The input file is not correctly formatted. Please provide a .txt with a column with SAMPLE_ID, POPULATION_ID, SUPERPOPULATION_ID, GENDER (Optional)')
            print('Exit...')
            sys.exit()
            return None, None, None, None, None, None, None, None
        #Add info of first line
        sample_to_pop[line[0]] = line[1]
        pop_to_superpop[line[1]] = line[2]
        
        superpop_to_pop[line[2]] = [line[1]]
        pop_to_sample[line[1]] = line[0]

        all_samples.add(line[0])
        all_pop.add(line[1])
        all_superpop.add(line[2])

        try:
            gender_sample[line[0]] = line[3]
        except:
            gender_sample[line[0]] = 'n/a'

        for line in in_file:
            line = line.strip().split('\t')
            sample_to_pop[line[0]] = line[1]
            pop_to_superpop[line[1]] = line[2]

            try:
                superpop_to_pop[line[2]].add(line[1])
            except :
                superpop_to_pop[line[2]] = set()
                superpop_to_pop[line[2]].add(line[1])
            try:
                pop_to_sample[line[1]].add(line[0])
            except:
                pop_to_sample[line[1]] = set()
                pop_to_sample[line[1]].add(line[0])
            
            all_samples.add(line[0])
            all_pop.add(line[1])
            all_superpop.add(line[2])

            try:
                gender_sample[line[0]] = line[3]
            except:
                gender_sample[line[0]] = 'n/a'
    return sample_to_pop, pop_to_superpop, superpop_to_pop, pop_to_sample, list(all_samples), list(all_pop), list(all_superpop), gender_sample

if __name__ == '__main__':
    loadSampleAssociation(sys.argv[1])

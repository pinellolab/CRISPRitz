#!/usr/bin/env python

#Dato in input un target.txt, la funzione calcola, per ogni diverso sample, il suo numero di occorrenze. Come prima riga abbiamo la guida + il numero di targets che hanno almeno un sample
#NON CANCELLARE
#Aggiunte colonne pam dis e pam creation
import os
import subprocess
import sys
import pandas as pd
# argv1  is result file, from top1 with expanded samples
# argv2 is job_id
# argv3 is type genome 'ref', 'var', 'both'
# argv4 is guides file
#NOTE Function only with vcf of HG38 and population info from the 20130606_sample_info.xlsx file

# NOTE PAM Disruptio  -> removed PAM Disruption calculation

pop_file = pd.read_excel(os.path.dirname(os.path.realpath(__file__)) + '/20130606_sample_info.xlsx')
all_samples = pop_file.Sample.to_list()
all_pop = pop_file.Population.to_list()
all_gender = pop_file.Gender.to_list()
dict_pop = dict()
gender_sample = dict()
for  pos, i in enumerate(all_samples):
    dict_pop[i] = all_pop[pos]
    gender_sample[i] = all_gender[pos]
population_1000gp = {'CHB':'EAS', 'JPT':'EAS', 'CHS':'EAS', 'CDX':'EAS', 'KHV':'EAS',
                    'CEU':'EUR', 'TSI':'EUR', 'FIN':'EUR', 'GBR':'EUR', 'IBS':'EUR',
                    'YRI':'AFR', 'LWK':'AFR', 'GWD':'AFR', 'MSL':'AFR', 'ESN':'AFR', 'ASW':'AFR', 'ACB':'AFR',
                    'MXL':'AMR', 'PUR':'AMR', 'CLM':'AMR', 'PEL':'AMR',
                    'GIH':'SAS', 'PJL':'SAS', 'BEB':'SAS', 'STU':'SAS', 'ITU':'SAS'
}
# Each guide has a dictionary, with samples as keys. Each sample (HG0096) has a list -> [Total targets, Var_uniq targets]
guides_dict = dict()
guides_dict_total = dict()  #contains total_sample_per_guide, did not merge the dict because didn't have time to do it
guides_population_targets = dict()
guides_superpopulation_targets = dict()
count_creation = dict() #{GUIDE1 -> {SAMPLE1 -> 0, SAMPLE2->1}, GUIDE2 -> {SAMPLE1 -> 0, SAMPLE2->1}}
current_chr_pos = '0'
with open(sys.argv[1]) as sample_file: #, open(sys.argv[3] + '.summary_by_samples.' + guide + '.txt', 'w+') as result:
    for line in sample_file:
        if '#' in line:
            continue
        line = line.strip().split('\t')
        if 'n' != line[-3]:
            guide = line[1].replace('-','')
            if guide not in guides_dict:
                guides_dict[guide] = dict()
                guides_dict_total[guide] = 0
                guides_population_targets[guide] = dict()       #{GUIDE -> {POP->count,POP2 -> count}, GUIDE2 -> {POP1 -> count}}
                guides_superpopulation_targets[guide] = dict()  #{GUIDE -> {SUPOP->count,SUPOP2 -> count}, GUIDE2 -> {SUPOP1 -> count}}
                count_creation[guide] = dict()
            words = line[-3].split(',')
            
            if current_chr_pos != guide + line[3]+line[4]:          
                guides_dict_total[guide] += 1                        #This value is skipped on app_v6
                current_chr_pos = guide + line[3]+line[4]
                checked_pop = []
                checked_superpop = []
                        
            for word in words:
                try:
                    guides_dict[guide][word][0] += 1
                except:
                    guides_dict[guide][word] = [1,0]
                if word not in count_creation[guide]:
                    count_creation[guide][word] = [0,0]
                
                if line[10] != 'n':
                    count_creation[guide][word][0] +=1

                if sys.argv[3] == 'both':
                    guides_dict[guide][word][1] += 1    #No more distinction for unique, add +1 to the [1] position meaning that the target was found in enriched
                    
                    if line[10] != 'n':
                        count_creation[guide][word][1] +=1

                #NOTE se voglio contare solo gli y il codice sotto lo idento una volta
                if dict_pop[word] not in checked_pop:       #I wanna count the target only once even if multiple samples of the same pop are in line[-3]
                    checked_pop.append(dict_pop[word])
                    try:
                        guides_population_targets[guide][dict_pop[word]] += 1
                    except:     #First time seeing this population
                        guides_population_targets[guide][dict_pop[word]] = 1  
                if  population_1000gp[dict_pop[word]] not in checked_superpop:
                    checked_superpop.append(population_1000gp[dict_pop[word]])
                    try:
                        guides_superpopulation_targets[guide][population_1000gp[dict_pop[word]]] += 1
                    except:
                        guides_superpopulation_targets[guide][population_1000gp[dict_pop[word]]] = 1

with open(sys.argv[4], 'r') as g_file:
    for line in g_file:
        line = line.strip()
        if line not in guides_dict:
            with open(sys.argv[2] + '.summary_by_samples.' + line + '.txt', 'w+') as result:
                result.write('No samples found with ' + line + ' guide')

#Add class from jobid.SampleClasses.txt
total_classes = dict()
dict_sample_classes = dict()
try:
    with open (sys.argv[2] + '.SampleClasses.txt') as sample_classes:
        headerclass = next(sample_classes).strip().split('\t')[1:]      #List of guides
        for g in headerclass:
            dict_sample_classes[g] = dict()
            total_classes[g] = dict()
        for line in sample_classes:
            if 'Total for' in line:
                line = line.strip().split('\t')
                for g, el in enumerate(line[1:]): 
                    total_classes[headerclass[g]] = el
                break
            line = line.strip().split('\t')
            for g, el in enumerate(line[1:]):             #iterate over the row, each cell is a class value for the specified guide (column)
                dict_sample_classes[headerclass[g]][line[0]] = el               #line[0] is sample
except FileNotFoundError as e:
    pass        #File SampleCLass does not exists

#Adjust general_target_count
res_tmp = []
try:
    with open(sys.argv[2] + '.general_target_count.txt') as gtc:
        header_tmp = next(gtc).strip().split('\t')
        header_tmp[1] = 'On-Targets Reference'
        header_tmp.insert(2,'Samples in Class 0 - 0+ - 1 - 1+')
        res_tmp.append('\t'.join(header_tmp))
        for line in gtc:
            line = line.strip().split('\t')
            line[1] = line[1].split('(')[-1].split('-')[0]    #Save only reference count
            line.insert(2, total_classes[line[0]])
            res_tmp.append('\t'.join(line))
    with open(sys.argv[2] + '.general_target_count.txt', 'w+') as gtc:
        gtc.write('\n'.join(res_tmp))
except FileNotFoundError as e:
    pass        #NO general_target_count_file


for k in guides_dict.keys():
    with open(sys.argv[2] + '.summary_by_samples.' + k + '.txt', 'w+') as result:
        # result.write(k + '\t' + str(guides_dict_total[k]) + '\n')       #this line is not used in the dash application
        if sys.argv[3] == 'both':
            result.write('\t'.join(['Sample', 'Gender', 'Population', 'Super Population',  'Targets in Reference', 'Targets in Enriched', 'Targets in Population', 'Targets in Super Population', 'PAM Creation', 'Class']) + '\n')
            for i in all_samples:
                if i not in guides_dict[k]:
                    guides_dict[k][i] = [0,0]
                    count_creation[k][i] = [0,0]
                if dict_pop[i] not in guides_population_targets[k]:
                    guides_population_targets[k][dict_pop[i]] = 0
                if population_1000gp[dict_pop[i]] not in guides_superpopulation_targets[k]:
                    guides_superpopulation_targets[k][population_1000gp[dict_pop[i]]] = 0


                result.write(i + '\t' + gender_sample[i] + '\t' + dict_pop[i] + '\t' + population_1000gp[dict_pop[i]] +'\t' + str(guides_dict[k][i][0]) + '\t' + str(guides_dict[k][i][1]) +
                            '\t' + str(guides_population_targets[k][dict_pop[i]]) + '\t' + 
                            str(guides_superpopulation_targets[k][population_1000gp[dict_pop[i]]])  + 
                            # '\t' + 
                            # str(count_disruption[k][i][0]) + 
                            '\t' + str(count_creation[k][i][1]) +
                            '\t' + str(dict_sample_classes[k][i]) + '\n')
        else:
            result.write('\t'.join(['Sample', 'Gender', 'Population', 'Super Population',  'Targets in Reference', 'Targets in Population', 'Targets in Super Population', 'PAM Creation']) + '\n')
            for i in all_samples:
                if i not in guides_dict[k]:
                    guides_dict[k][i] = [0,0]
                    count_creation[k][i] = [0,0]
                if dict_pop[i] not in guides_population_targets[k]:
                    guides_population_targets[k][dict_pop[i]] = 0
                if population_1000gp[dict_pop[i]] not in guides_superpopulation_targets[k]:
                    guides_superpopulation_targets[k][population_1000gp[dict_pop[i]]] = 0

                result.write(i + '\t' + gender_sample[i] + '\t' + dict_pop[i] + '\t' + population_1000gp[dict_pop[i]] +'\t' + str(guides_dict[k][i][0]) +
                            '\t' + str(guides_population_targets[k][dict_pop[i]]) + '\t' + str(guides_superpopulation_targets[k][population_1000gp[dict_pop[i]]])  + 
                            # '\t' + 
                            # str(count_disruption[k][i][0]) + 
                            '\t' + str(count_creation[k][i][1]) +'\n')


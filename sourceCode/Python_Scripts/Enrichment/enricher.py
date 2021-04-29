#!/usr/bin/env python

import sys
import re
import time
import textwrap
import os
import json
import pandas as pd
import gzip

start_time = time.time()

print("READING VCF FILE AND CHROMOSOME")

altFile = sys.argv[1]  # file with variations
genomeFile = sys.argv[2]  # genome file
dir_enr_name = sys.argv[3] + '_enriched' #name of directory for saving files
inAltFile = gzip.open(altFile, "rt")  # variations file open
inGenomeFile = open(genomeFile, "r")  # genome file open
vcfName = os.path.basename(sys.argv[5])

#read fasta header (chr_name)
genomeHeader = inGenomeFile.readline()
currentChr = genomeHeader[1:].strip()
if 'chr' not in currentChr:
    currentChr = 'chr'+currentChr

#read VCF header
VCFheader=inAltFile.readline().split('\t')

#all the fasta content is read and stored in a single line with \n removed and uppercased to avoid possible conflicts
genomeStr = inGenomeFile.read()
genomeStr = genomeStr.replace('\n', '')
genomeStr = genomeStr.upper()

#genome str is listed to be modified faster and easier
genomeList = list(genomeStr)

iupac_code = {}

iupac_code.setdefault('R', []).append('AG')
iupac_code.setdefault('R', []).append('GA')

iupac_code.setdefault('Y', []).append('CT')
iupac_code.setdefault('Y', []).append('TC')

iupac_code.setdefault('S', []).append('GC')
iupac_code.setdefault('S', []).append('CG')

iupac_code.setdefault('W', []).append('AT')
iupac_code.setdefault('W', []).append('TA')

iupac_code.setdefault('K', []).append('GT')
iupac_code.setdefault('K', []).append('TG')

iupac_code.setdefault('M', []).append('AC')
iupac_code.setdefault('M', []).append('CA')

iupac_code.setdefault('B', []).append('CGT')
iupac_code.setdefault('B', []).append('GCT')
iupac_code.setdefault('B', []).append('TGC')
iupac_code.setdefault('B', []).append('GTC')
iupac_code.setdefault('B', []).append('CTG')
iupac_code.setdefault('B', []).append('TCG')

iupac_code.setdefault('D', []).append('AGT')
iupac_code.setdefault('D', []).append('GAT')
iupac_code.setdefault('D', []).append('TAG')
iupac_code.setdefault('D', []).append('ATG')
iupac_code.setdefault('D', []).append('GTA')
iupac_code.setdefault('D', []).append('TGA')

iupac_code.setdefault('H', []).append('ACT')
iupac_code.setdefault('H', []).append('CAT')
iupac_code.setdefault('H', []).append('TCA')
iupac_code.setdefault('H', []).append('ATC')
iupac_code.setdefault('H', []).append('CTA')
iupac_code.setdefault('H', []).append('TAC')

iupac_code.setdefault('V', []).append('ACG')
iupac_code.setdefault('V', []).append('CAG')
iupac_code.setdefault('V', []).append('GAC')
iupac_code.setdefault('V', []).append('AGC')
iupac_code.setdefault('V', []).append('CGA')
iupac_code.setdefault('V', []).append('GCA')

iupac_code.setdefault('N', []).append('ACGT')
iupac_code.setdefault('N', []).append('CAGT')
iupac_code.setdefault('N', []).append('GACT')
iupac_code.setdefault('N', []).append('AGCT')
iupac_code.setdefault('N', []).append('CGAT')
iupac_code.setdefault('N', []).append('GCAT')
iupac_code.setdefault('N', []).append('GCTA')
iupac_code.setdefault('N', []).append('CGTA')
iupac_code.setdefault('N', []).append('TGCA')
iupac_code.setdefault('N', []).append('GTCA')
iupac_code.setdefault('N', []).append('CTGA')
iupac_code.setdefault('N', []).append('TCGA')
iupac_code.setdefault('N', []).append('TAGC')
iupac_code.setdefault('N', []).append('ATGC')
iupac_code.setdefault('N', []).append('GTAC')
iupac_code.setdefault('N', []).append('TGAC')
iupac_code.setdefault('N', []).append('AGTC')
iupac_code.setdefault('N', []).append('GATC')
iupac_code.setdefault('N', []).append('CATG')
iupac_code.setdefault('N', []).append('ACTG')
iupac_code.setdefault('N', []).append('TCAG')
iupac_code.setdefault('N', []).append('CTAG')
iupac_code.setdefault('N', []).append('ATCG')
iupac_code.setdefault('N', []).append('TACG')

iupac_code_scomposition = {
          "R":"AG",
          "Y":"CT",
          "S":"GC",
          "W":"AT",
          "K":"GT",
          "M":"AC",
          "B":"CGT",
          "D":"AGT",
          "H":"ACT",
          "V":"ACG",
          "r":"AG",
          "y":"CT",
          "s":"GC",
          "w":"AT",
          "k":"GT",
          "m":"AC",
          "b":"CGT",
          "d":"AGT",
          "h":"ACT",
          "v":"ACG",
          "A":"A",
          "T":"T",
          "C":"C",
          "G":"G",
          "a":"A",
          "t":"T",
          "c":"C",
          "g":"G",
          'N':'ATGC'
        }


def SNPsProcess(line):
    if len(line[3]) > 1: #se len ref > 1 significa che non è snp
        return
    line[1] = str(int(line[1])-1)  # taaac per corregere conta posizione da zero/uno
    altAlleles = line[4].strip().split(',') #list di tutti gli alleli alternative
    referenceNucleotide = genomeStr[int(line[1])] #nucleotide reference preso dalla stringa del fasta
    iupac_value = set() #set che contiene nt reference + alt snps
    iupac_value.add(referenceNucleotide)
    for alt in altAlleles: #ciclo sugli alt nt
        if len(alt) > 1: #se lun di alt nt è > 1 skippo perchè non è snp
            continue
        else:
            iupac_value.add(alt) #add del nt alt al set di iupac
    iupac_value = ''.join(iupac_value) #converto il set in una string
    for key, value in iupac_code.items(): #cerco iupac value nel dict per avere il corrispettivo IUPAC nt
        if iupac_value in value:
            genomeList[int(line[1])] = str(key) #inserisco il IUPAC nt nel genoma enriched


# def SNPsProcess(line):
#     #x = line.strip().split('\t') #split alt file to add snps to genome
#     x = line
#     #del x[1] #remove the unnecessary rsID from the split list
#     x[1] = str(int(x[1])-1)  # taaac
#     if (',' in x[4]) and (len(x[3]) == 1) and ('>' not in x[4]) and (len(x[4]) < 6):
#         altstr = str(x[4])
#         k = altstr.split(',')
#         if (len(x[4]) == 3):
#             if (str(genomeList[int(x[1])])) not in k:
#                 original = genomeList[int(x[1])]
#             else:
#                 original = ''
#             snp = k[0]+k[1]
#             iupacvalue = str(original+snp)  #TODO if original is a single nucleotide, use iupac_code_scomposition[original]
#             for key, value in iupac_code.items():
#                 if iupacvalue in value:
#                     genomeList[int(x[1])] = str(key)
#         elif (len(x[4]) == 5):
#             if (len(k[0]) == 1) and (len(k[1]) == 1) and (len(k[2]) == 1):
#                 if (str(genomeList[int(x[1])])) not in k:
#                     original = genomeList[int(x[1])]
#                 else:
#                     original = ''
#                 original = genomeList[int(x[1])]
#                 snp = k[0]+k[1]+k[2]
#                 iupacvalue = str(original+snp)  #TODO if original is a single nucleotide, use iupac_code_scomposition[original]
#                 for key, value in iupac_code.items():
#                     if iupacvalue in value:
#                         genomeList[int(x[1])] = str(key)
#     elif (',' not in x[4]) and (len(x[3]) == 1) and ('>' not in x[4]) and (len(x[4]) == 1):
#         original = iupac_code_scomposition[genomeList[int(x[1])]]
#         snp = x[4]
#         iupacvalue = str(original+snp)
#         for key, value in iupac_code.items():
#             if iupacvalue in value:
#                 genomeList[int(x[1])] = str(key)

def chromosomeSave():
    genomeStr = "".join(genomeList)
    os.chdir("./SNPs_genome/")
    if not (os.path.isdir(dir_enr_name)):
        os.mkdir(dir_enr_name)
    outFile = open(dir_enr_name + '/' + genomeHeader[1:(len(genomeHeader)-1)]+'.enriched'+'.fa', 'w')
    outFile.write(genomeHeader+genomeStr+'\n')
    outFile.close()

def add_to_dict_snps(line, pos_AF):
    list_samples = []
    list_chars = []
    if len(line[3]) == 1 and len(line[4]) == 1:
        for pos, i in enumerate(line[9:]):          #if sample has 1|1 0|1 or 1|0, #NOTE may change for different vcf
            if ('1' in i.split(':')[0]):
                list_samples.append(VCFheader[ pos + 9]+'('+i.split(':')[0]+')')
        
        # print(list_samples)
        
        chr_pos_string = currentChr + ',' + line[1] #chr,position
        #Add in last two position the ref and alt nucleotide, eg: chrX,100 -> sample1,sample5,sample10;A,T;rsID100;0.01
        #If no sample was found, the dict is chrX,100 -> ;A,T;rsID100;0.01
        rsID = line[2]
        list_chars.append(line[3])
        list_chars.append(line[4])
        af = line[7].split(";")[pos_AF][3:] 
        if len(list_samples) > 0:
            chr_dict_snps[chr_pos_string] = ','.join(sorted(list_samples)) + ';' + ','.join(list_chars) + ";" + rsID + ";" + af
        else:
            chr_dict_snps[chr_pos_string] = ';' + ','.join(list_chars) + ";" + rsID + ";" + af #None
        print(chr_dict_snps[chr_pos_string])
    elif len(line[3]) == 1:
        variants = line[4].split(",")
        snps = []
        values_for_allele_info = []
        for value, var in enumerate(variants):
            if len(var) == 1: 
                snps.append(var)
                values_for_allele_info.append(value+1) #save number corresponding to snp (1,2,3...)
        dict_of_lists_samples = {}
        for snp in snps:
            dict_of_lists_samples[snp] = []
        if len(snps) > 0:		
            for pos, sample in enumerate(line[9:]):
                for idx, value in enumerate(values_for_allele_info):
                    if str(value) in sample.split(':')[0]:
                        dict_of_lists_samples[snps[idx]].append(VCFheader[ pos + 9]+'('+sample.split(':')[0]+')') #add to correct entry of dict the sample with such snp
                        break
                    
            chr_pos_string = currentChr + ',' + line[1]
            rsID = line[2].split(',')
            af = line[7].split(";")[pos_AF][3:].split(',')

            final_entry = []
            for idx, snp in enumerate(snps):
                list_chars = [line[3]]
                list_chars.append(snp)
                if len(dict_of_lists_samples[snp]) > 0:
                    final_entry.append(','.join(sorted(dict_of_lists_samples[snp])) + ';' + ','.join(list_chars) + ";" + rsID[0] + ";" + af[values_for_allele_info[idx]-1])
                else:
                    final_entry.append(';' + ','.join(list_chars) + ";" + rsID[0] + ";" + af[values_for_allele_info[idx]-1])
            chr_dict_snps[chr_pos_string] = '/'.join(final_entry)
        print('/'.join(final_entry))

def dictSave():
    #os.chdir("./SNPs_genome/")
    with open('my_dict_' + currentChr + '.json', 'w') as f:
        json.dump(chr_dict_snps, f) 


def indel_to_fasta(line, id_indel, pos_AF, start_fake_pos):
    if (len(line[3]) > 1 or len(line[4]) > 1) and '<' not in line[3]:
        #print(''.join(line[3:5]))
        indels = []
        values_for_allele_info = []
        if ',' in line[4]:
            splitted = line[4].split(',')
            for value, s in enumerate(splitted):
                if (len(s) != len(line[3]) or (len(s) > 1 and len(line[3]) > 1)) and '<' not in s:
                    indels.append(s)
                    values_for_allele_info.append(value+1)
        elif (len(line[4]) != len(line[3]) or (len(line[4]) > 1 and len(line[3]) > 1)) and '<' not in line[4]:
            indels.append(line[4])
            values_for_allele_info.append(1)
        #print(search_sample_value, line[3], line[4])
        if len(indels) > 0:

            dict_of_lists_samples = {}
            for indel in indels:
                dict_of_lists_samples[indel] = []

            for pos, sample in enumerate(line[9:]):          #if sample has 1|1 0|1 or 1|0, #NOTE may change for different vcf
                for idx, value in enumerate(values_for_allele_info):
                    if str(value) in sample.split(':')[0]:
                        dict_of_lists_samples[indels[idx]].append(VCFheader[ pos + 9]) #add to correct entry of dict the sample with such snp
                        break
            if len(indels) > 0:
                rsID = line[2].split(',')
                af = line[7].split(";")[pos_AF][3:].split(',')

                start_position = int(line[1])-26
                end_position = int(line[1])+26+len(line[3])
                sub_fasta = genomeStr[start_position:end_position]
                for idx, indel in enumerate(indels):
                    if len(dict_of_lists_samples[indel]) > 0:
                        indel_info = f"{currentChr}_{line[1]}_{line[3]}_{indel}"
                        
                        #sub_fasta[25] = line[3] 
                        #sub_fasta = ''.join(sub_fasta)
                        sub_fasta = sub_fasta[0:25] + re.sub(line[3], indel, sub_fasta[25:], 1, flags=re.IGNORECASE) 
                        
                        #fasta_out.write(f'>{currentChr}_{start_position}-{end_position}_{id_indel}\n')
                        list_fasta_indels.append(sub_fasta+'\n'+"N\n")
                        
                        refseq = genomeStr[start_position:start_position+len(sub_fasta)]
                        end_fake_pos = start_fake_pos + len(sub_fasta)#(end_position - start_position)
                        
                        log_indels.append([f"{currentChr}_{start_position}-{end_position}_{id_indel}", ",".join(dict_of_lists_samples[indel]), rsID[0], af[values_for_allele_info[idx]-1], indel_info, f"{start_fake_pos},{end_fake_pos}", refseq])
                    
                        id_indel += 1
                        start_fake_pos = end_fake_pos + 1 

    return id_indel, start_fake_pos

def logIndelsSave():
    fasta_out.write(''.join(list_fasta_indels))
    fasta_out.close()
    df = pd.DataFrame(log_indels, columns=['CHR','SAMPLES','rsID','AF','indel', 'FAKEPOS', 'refseq'])
    df.to_csv('log' + currentChr + '.txt', index=False, sep='\t')
    #dict_indels = df.to_dict(orient='index')
    #with open('log' + currentChr + '.json', 'w') as f:
    #    json.dump(log_indels, f) 

print('START ENRICHMENT WITH SNVs AND SVs')

if sys.argv[4] == 'yes':
    chr_dict_snps = dict()
    log_indels = []
    id_indel = 1
    start_fake_pos = 0
    if not os.path.isdir(f"fake_{vcfName}_{currentChr}"):
        os.mkdir(f"fake_{vcfName}_{currentChr}")
    output_dir_indels = os.path.abspath(f"fake_{vcfName}_{currentChr}")
    fasta_out = open(f'{output_dir_indels}/fake{currentChr}.fa', 'w+')
    fasta_out.write(f">fake{currentChr}\n")
    list_fasta_indels = []

first_line = True
for line in inAltFile:
    if ('#CHROM') in line:
        VCFheader = line.strip().split('\t')   #Save this header for retrieving sample id
        break
for line in inAltFile:
    line = line.strip().split('\t')
    if first_line:
        first_line = False
        splitted = line[7].split(";")
        for pos, ele in enumerate(splitted):
            if ele[0:2] == "AF":
                pos_AF = pos
                break
    if line[6] != 'PASS':
        continue
    if sys.argv[4] == 'yes': #if true do all the creation, if false do only genome enrichment with SNPs
        #print(sys.argv[4])
        add_to_dict_snps(line, pos_AF)
        id_indel, start_fake_pos = indel_to_fasta(line, id_indel, pos_AF, start_fake_pos)
    SNPsProcess(line)
    

#saving chr after enrichment with snps
chromosomeSave()
if sys.argv[4] == 'yes':
    dictSave()
    logIndelsSave()

#REMOVED INDELS CREATION AND SUBSTITUTED WITH NEW SCRIPT

#step to generate fasta files with indels and snps, now INDELS are treated separately, this works only for single sample VCF
# for line in inAltFile:
#     # if 'PASS' not in line: #skip line with QUAL inferior to PASS
#     #     continue
#     x = line.strip().split('\t')
#     x[0] = str(int(x[0])-1)  # taaac
#     if (',' not in x[2]) and (',' not in x[1]) and ('>' not in x[2]) and (len(x[1]) == 1) and (len(x[2]) > 1):
#         genomeList[int(x[0])] = str(x[2])
#     elif (',' not in x[2]) and (',' not in x[1]) and ('>' not in x[2]) and (len(x[1]) > 1) and (len(x[2]) == 1):
#         point = int(x[0])
#         fine = point+len(x[1])
#         genomeList[point] = str(x[2])
#         genomeList[point+1:fine] = [''] * (len(x[1])-1)

# genomeStr = "".join(genomeList)

# os.chdir("../")
# os.chdir("./INDELs_genome/")
# if not (os.path.isdir(dir_enr_name)):
#     os.mkdir(dir_enr_name)
# outfile = open(dir_enr_name + '/' + genomeHeader[1:(len(genomeHeader)-1)]+'.indels'+'.fa', 'w')
# outfile.write(genomeHeader+genomeStr+'\n')

print("DONE IN ",int(time.time()-start_time),' seconds')

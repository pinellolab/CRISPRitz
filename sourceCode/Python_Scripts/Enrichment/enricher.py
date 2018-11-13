#!/usr/bin/python3

import sys
import re
import time
import textwrap


start_time = time.time()

##print "READING"

altFile = sys.argv[1]                 #file with variations
genomeFile = sys.argv[2]               #genome file         

inAltFile = open(altFile, "r").readlines()          #variations file open
inGenomeFile = open(genomeFile, "r")      #genome file open                  
      #outfile open                         

##print('start lettura e conversioni')

genomeHeader = inGenomeFile.readline()

genomeStr = inGenomeFile.read()
#genomeStr = genomeStr.replace('\n','')
genomeStr = genomeStr.upper()

genomeList = list(genomeStr)



iupac_code = {}

iupac_code.setdefault('R',[]).append('AG')
iupac_code.setdefault('R',[]).append('GA')

iupac_code.setdefault('Y',[]).append('CT')
iupac_code.setdefault('Y',[]).append('TC')

iupac_code.setdefault('S',[]).append('GC')
iupac_code.setdefault('S',[]).append('CG')

iupac_code.setdefault('W',[]).append('AT')
iupac_code.setdefault('W',[]).append('TA')

iupac_code.setdefault('K',[]).append('GT')
iupac_code.setdefault('k',[]).append('TG')

iupac_code.setdefault('M',[]).append('AC')
iupac_code.setdefault('M',[]).append('CA')

iupac_code.setdefault('B',[]).append('CGT')
iupac_code.setdefault('B',[]).append('GCT')
iupac_code.setdefault('B',[]).append('TGC')
iupac_code.setdefault('B',[]).append('GTC')
iupac_code.setdefault('B',[]).append('CTG')
iupac_code.setdefault('B',[]).append('TCG')

iupac_code.setdefault('D',[]).append('AGT')
iupac_code.setdefault('D',[]).append('GAT')
iupac_code.setdefault('D',[]).append('TAG')
iupac_code.setdefault('D',[]).append('ATG')
iupac_code.setdefault('D',[]).append('GTA')
iupac_code.setdefault('D',[]).append('TGA')

iupac_code.setdefault('H',[]).append('ACT')
iupac_code.setdefault('H',[]).append('CAT')
iupac_code.setdefault('H',[]).append('TCA')
iupac_code.setdefault('H',[]).append('ATC')
iupac_code.setdefault('H',[]).append('CTA')
iupac_code.setdefault('H',[]).append('TAC')

iupac_code.setdefault('V',[]).append('ACG')
iupac_code.setdefault('V',[]).append('CAG')
iupac_code.setdefault('V',[]).append('GAC')
iupac_code.setdefault('V',[]).append('AGC')
iupac_code.setdefault('V',[]).append('CGA')
iupac_code.setdefault('V',[]).append('GCA')

iupac_code.setdefault('N',[]).append('ACGT')
iupac_code.setdefault('N',[]).append('CAGT')
iupac_code.setdefault('N',[]).append('GACT')
iupac_code.setdefault('N',[]).append('AGCT')
iupac_code.setdefault('N',[]).append('CGAT')
iupac_code.setdefault('N',[]).append('GCAT')
iupac_code.setdefault('N',[]).append('GCTA')
iupac_code.setdefault('N',[]).append('CGTA')
iupac_code.setdefault('N',[]).append('TGCA')
iupac_code.setdefault('N',[]).append('GTCA')
iupac_code.setdefault('N',[]).append('CTGA')
iupac_code.setdefault('N',[]).append('TCGA')
iupac_code.setdefault('N',[]).append('TAGC')
iupac_code.setdefault('N',[]).append('ATGC')
iupac_code.setdefault('N',[]).append('GTAC')
iupac_code.setdefault('N',[]).append('TGAC')
iupac_code.setdefault('N',[]).append('AGTC')
iupac_code.setdefault('N',[]).append('GATC')
iupac_code.setdefault('N',[]).append('CATG')
iupac_code.setdefault('N',[]).append('ACTG')
iupac_code.setdefault('N',[]).append('TCAG')
iupac_code.setdefault('N',[]).append('CTAG')
iupac_code.setdefault('N',[]).append('ATCG')
iupac_code.setdefault('N',[]).append('TACG')

#print('start ciclo')


for line in inAltFile:
    x = line.split(' ')
    if (',' in x[2]) and (len(x[1])==1) and ('>' not in x[2]) and (len(x[2])<6) and ('\n' not in genomeList[int(x[0])]):
        altstr = str(x[2])
        k = altstr.split(',')
        if (len(x[2])==3):
            if (str(genomeList[int(x[0])])) not in k:
                original = genomeList[int(x[0])]
            else:
                original = ''
            snp = k[0]+k[1]
            iupacvalue = str(original+snp)
            for key, value in iupac_code.items():
                if iupacvalue in value:
                    genomeList[int(x[0])] = str(key)
        elif (len(x[2])==5):
            if (len(k[0])==1) and (len(k[1])==1) and (len(k[2])==1):
                if (str(genomeList[int(x[0])])) not in k:
                    original = genomeList[int(x[0])]
                else:
                    original = ''
                original = genomeList[int(x[0])]
                snp = k[0]+k[1]+k[2]
                iupacvalue = str(original+snp)
            for key, value in iupac_code.items():
                if iupacvalue in value:
                    genomeList[int(x[0])] = str(key)
    elif (',' not in x[2]) and (len(x[1])==1) and ('>' not in x[2]) and (len(x[2])==1) and ('\n' not in genomeList[int(x[0])]):
        original = genomeList[int(x[0])]
        snp = x[2]
        iupacvalue = str(original+snp)
        for key, value in iupac_code.items():
            if iupacvalue in value:
                genomeList[int(x[0])] = str(key)


genomeStr = "".join(genomeList)

outFile = open(genomeHeader[1:(len(genomeHeader)-1)]+'.enriched'+'.fa', 'w')
outFile.write(genomeHeader+genomeStr)

# genomeIndel = []
# start=0
# end=0

for line in inAltFile:
    x = line.split(' ')
    if (',' not in x[2]) and (',' not in x[1]) and ('>' not in x[2]) and (len(x[1])==1) and (len(x[2])>1) and ('\n' not in genomeList[int(x[0]):(int(x[0])+len(x[2]))]):
        #print ("fatto 1")
        # genomeList[int(x[0])] = str(x[2])
        # end=(int(x[0])-1)
        # genomeIndel.append(genomeStr[start:end] + str(x[2]))
        # start=(int(x[0])+1
        #print(genomeList[int(x[0])])
        genomeList[int(x[0])] = str(x[2])
        #print(genomeList[int(x[0])])
    elif (',' not in x[2]) and (',' not in x[1]) and ('>' not in x[2]) and (len(x[1])>1) and (len(x[2])==1) and ('\n' not in genomeList[int(x[0]):(int(x[0])+len(x[1]))]):
        #print ("fatto 2")
        point = int(x[0])
        fine = point+len(x[1])
        #genomeList[point:(point+sublen)]=str(x[2])
        # end=(int(x[0])-1)
        # genomeIndel.append(genomeStr[start:end] + str(x[2]))
        # start=(int(x[0])+len(x[1]))+1
        #print(genomeList[point:fine])
        #print(genomeList[point])
        genomeList[point] = str(x[2])
        #print(genomeList[point])
        #print(genomeList[point+1:fine])
        genomeList[point+1:fine] = [""] * (len(x[1])-1)
        #print(genomeList[point+1:fine])

    
genomeStr = "".join(genomeList)

outfile = open(genomeHeader[1:(len(genomeHeader)-1)]+'.indels'+'.fa', 'w')
outfile.write(genomeHeader+genomeStr)

print("Done")


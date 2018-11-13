#!/usr/bin/python3

from intervaltree import Interval, IntervalTree

import sys
import time
import concurrent.futures


print ("READING")

#exonsFile = sys.argv[1]                 #exonsfile path
#intronsFile = sys.argv[2]               #intronsfile path
guidesFile = sys.argv[6]         
resultFile = sys.argv[7]                #resultsfile path                                  
annotatedFile = sys.argv[8]             #annotatedfile path
#inExons = open(exonsFile, "r")          #exonsfile open
#inIntrons = open(intronsFile, "r")      #intronsfile open                  
inResult = open(resultFile, "r")        #resultfile open
outFile = open(annotatedFile + '.txt', 'w')      #outfile open  
inGuides = open(guidesFile, "r")   


start_time = time.time()


print ("CREATING TREES")
exonsFile = sys.argv[1]
intronsFile = sys.argv[2]
promotersFile = sys.argv[3]
ctcfFile = sys.argv[4]
dnaseFile = sys.argv[5]

d_exons_0MM = {}
d_exons_1MM = {}
d_exons_2MM = {}
d_exons_3MM = {}
d_exons_4MM = {}
d_exons_5MM = {}
d_exons_6MM = {}
d_exons_7MM = {}
d_exons_8MM = {}
d_exons_9MM = {}
d_exons_10MM = {}

d_introns_0MM = {}
d_introns_1MM = {}
d_introns_2MM = {}
d_introns_3MM = {}
d_introns_4MM = {}
d_introns_5MM = {}
d_introns_6MM = {}
d_introns_7MM = {}
d_introns_8MM = {}
d_introns_9MM = {}
d_introns_10MM = {}

d_promoters_0MM = {}
d_promoters_1MM = {}
d_promoters_2MM = {}
d_promoters_3MM = {}
d_promoters_4MM = {}
d_promoters_5MM = {}
d_promoters_6MM = {}
d_promoters_7MM = {}
d_promoters_8MM = {}
d_promoters_9MM = {}
d_promoters_10MM = {}

d_ctcf_0MM = {}
d_ctcf_1MM = {}
d_ctcf_2MM = {}
d_ctcf_3MM = {}
d_ctcf_4MM = {}
d_ctcf_5MM = {}
d_ctcf_6MM = {}
d_ctcf_7MM = {}
d_ctcf_8MM = {}
d_ctcf_9MM = {}
d_ctcf_10MM = {}

d_dnase_0MM = {}
d_dnase_1MM = {}
d_dnase_2MM = {}
d_dnase_3MM = {}
d_dnase_4MM = {}
d_dnase_5MM = {}
d_dnase_6MM = {}
d_dnase_7MM = {}
d_dnase_8MM = {}
d_dnase_9MM = {}
d_dnase_10MM = {}

exonsTree = IntervalTree()
intronsTree = IntervalTree()
promotersTree = IntervalTree()
ctcfTree = IntervalTree()
dnaseTree = IntervalTree()

for line in inGuides:
    x = line.split('\n')
    d_exons_0MM[str(x[0])] = 0
    d_exons_1MM[str(x[0])] = 0
    d_exons_2MM[str(x[0])] = 0
    d_exons_3MM[str(x[0])] = 0
    d_exons_4MM[str(x[0])] = 0
    d_exons_5MM[str(x[0])] = 0
    d_exons_6MM[str(x[0])] = 0
    d_exons_7MM[str(x[0])] = 0
    d_exons_8MM[str(x[0])] = 0
    d_exons_9MM[str(x[0])] = 0
    d_exons_10MM[str(x[0])] = 0

    d_introns_0MM[str(x[0])] = 0
    d_introns_1MM[str(x[0])] = 0
    d_introns_2MM[str(x[0])] = 0
    d_introns_3MM[str(x[0])] = 0
    d_introns_4MM[str(x[0])] = 0
    d_introns_5MM[str(x[0])] = 0
    d_introns_6MM[str(x[0])] = 0
    d_introns_7MM[str(x[0])] = 0
    d_introns_8MM[str(x[0])] = 0
    d_introns_9MM[str(x[0])] = 0
    d_introns_10MM[str(x[0])] = 0

    d_dnase_0MM[str(x[0])] = 0
    d_dnase_1MM[str(x[0])] = 0
    d_dnase_2MM[str(x[0])] = 0
    d_dnase_3MM[str(x[0])] = 0
    d_dnase_4MM[str(x[0])] = 0
    d_dnase_5MM[str(x[0])] = 0
    d_dnase_6MM[str(x[0])] = 0
    d_dnase_7MM[str(x[0])] = 0
    d_dnase_8MM[str(x[0])] = 0
    d_dnase_9MM[str(x[0])] = 0
    d_dnase_10MM[str(x[0])] = 0

    d_promoters_0MM[str(x[0])] = 0
    d_promoters_1MM[str(x[0])] = 0
    d_promoters_2MM[str(x[0])] = 0
    d_promoters_3MM[str(x[0])] = 0
    d_promoters_4MM[str(x[0])] = 0
    d_promoters_5MM[str(x[0])] = 0
    d_promoters_6MM[str(x[0])] = 0
    d_promoters_7MM[str(x[0])] = 0
    d_promoters_8MM[str(x[0])] = 0
    d_promoters_9MM[str(x[0])] = 0
    d_promoters_10MM[str(x[0])] = 0

    d_ctcf_0MM[str(x[0])] = 0
    d_ctcf_1MM[str(x[0])] = 0
    d_ctcf_2MM[str(x[0])] = 0
    d_ctcf_3MM[str(x[0])] = 0
    d_ctcf_4MM[str(x[0])] = 0
    d_ctcf_5MM[str(x[0])] = 0
    d_ctcf_6MM[str(x[0])] = 0
    d_ctcf_7MM[str(x[0])] = 0
    d_ctcf_8MM[str(x[0])] = 0
    d_ctcf_9MM[str(x[0])] = 0
    d_ctcf_10MM[str(x[0])] = 0


if(exonsFile != "no"):
    inExons = open(exonsFile, "r")
    for line in inExons:
        x = line.split('\t')                        
        exonsTree[int(x[1]):int(x[2])] = str(x[0])

if(intronsFile != "no"):
    inIntrons = open(intronsFile, "r") 
    for line in inIntrons:
        x = line.split('\t')
        intronsTree[int(x[1]):int(x[2])] = str(x[0])

if(promotersFile != "no"):
    inPromoters = open(promotersFile, "r") 
    for line in inPromoters:
        x = line.split('\t')
        promotersTree[int(x[1]):int(x[2])] = str(x[0])

if(ctcfFile != "no"):
    inCtcf = open(ctcfFile, "r") 
    for line in inCtcf:
        x = line.split('\t')
        ctcfTree[int(x[1]):int(x[2])] = str(x[0])

if(dnaseFile != "no"):
    inDnase = open(dnaseFile, "r") 
    for line in inDnase:
        x = line.split('\t')
        dnaseTree[int(x[1]):int(x[2])] = str(x[0])

print("TREES CREATION TIME: %s seconds" % (time.time() - start_time))

# print ("CREATING TREES")

# exons = IntervalTree()
# introns = IntervalTree()

# for line in inExons:
#     x = line.split('\t')                        #split the line
#     exons[int(x[1]):int(x[2])] = str(x[0])      #exons tree creation

# for line in inIntrons:
#     x = line.split('\t')                        #split the line
#     introns[int(x[1]):int(x[2])] = str(x[0])    #introns tree creation



start_time = time.time()

print ("ANNOTATING")


for line in inResult:
    x = line.split('\t')
    foundExons = sorted(exonsTree[int(x[2]):(int(x[2])+int(len(x[0]))+1)])              #exons interval query
    foundIntrons = sorted(intronsTree[int(x[2]):(int(x[2])+int(len(x[0]))+1)])          #introns interval query
    foundPromoters = sorted(promotersTree[int(x[2]):(int(x[2])+int(len(x[0]))+1)])
    foundCtcf = sorted(ctcfTree[int(x[2]):(int(x[2])+int(len(x[0]))+1)])
    foundDnase = sorted(dnaseTree[int(x[2]):(int(x[2])+int(len(x[0]))+1)])

    for item in range(0,len(foundExons)):
        if foundExons[item].data == str(x[1]):
            outFile.write(str(x[3]) + '\t' + str(foundExons[item].data) + '\t' + str(x[2]) + '\t' + str(int(x[2])+int(len(x[0]))) + '\t' + str(foundExons[item].begin) + '\t' 
            + str(foundExons[item].end) + '\t' + str(x[4]) + '\t' + "exon"+ '\t' +str(x[5]))       #file write of exons results
            if str(x[5]) == '0\n':
	            d_exons_0MM[str(x[0])] +=1
            elif str(x[5]) == '1\n':
	            d_exons_1MM[str(x[0])] +=1
            elif str(x[5]) == '2\n':
                d_exons_2MM[str(x[0])] +=1
            elif str(x[5]) == '3\n':
                d_exons_3MM[str(x[0])] +=1
            elif str(x[5]) == '4\n':
                d_exons_4MM[str(x[0])] +=1
            elif str(x[5]) == '5\n':
                d_exons_5MM[str(x[0])] +=1
            elif str(x[5]) == '6\n':
                d_exons_6MM[str(x[0])] +=1
            elif str(x[5]) == '7\n':
                d_exons_7MM[str(x[0])] +=1
            elif str(x[5]) == '8\n':
                d_exons_8MM[str(x[0])] +=1
            elif str(x[5]) == '9\n':
                d_exons_9MM[str(x[0])] +=1
            elif str(x[5]) == '10\n':
                d_exons_10MM[str(x[0])] +=1
    
    for item in range(0,len(foundIntrons)):
        if foundIntrons[item].data == str(x[1]):
            outFile.write(str(x[3]) + '\t' + str(foundIntrons[item].data) + '\t' + str(x[2]) + '\t' + str(int(x[2])+int(len(x[0]))) + '\t' + str(foundIntrons[item].begin) + '\t' 
            + str(foundIntrons[item].end) + '\t' + str(x[4]) + '\t' + "intron"+ '\t' +str(x[5]))  #file write of introns results
            if str(x[5]) == '0\n':
	            d_introns_0MM[str(x[0])] +=1
            elif str(x[5]) == '1\n':
	            d_introns_1MM[str(x[0])] +=1
            elif str(x[5]) == '2\n':
                d_introns_2MM[str(x[0])] +=1
            elif str(x[5]) == '3\n':
                d_introns_3MM[str(x[0])] +=1
            elif str(x[5]) == '4\n':
                d_introns_4MM[str(x[0])] +=1
            elif str(x[5]) == '5\n':
                d_introns_5MM[str(x[0])] +=1
            elif str(x[5]) == '6\n':
                d_introns_6MM[str(x[0])] +=1
            elif str(x[5]) == '7\n':
                d_introns_7MM[str(x[0])] +=1
            elif str(x[5]) == '8\n':
                d_introns_8MM[str(x[0])] +=1
            elif str(x[5]) == '9\n':
                d_introns_9MM[str(x[0])] +=1
            elif str(x[5]) == '10\n':
                d_introns_10MM[str(x[0])] +=1
            
    for item in range(0,len(foundPromoters)):
        if foundPromoters[item].data == str(x[1]):
            outFile.write(str(x[3]) + '\t' + str(foundPromoters[item].data) + '\t' + str(x[2]) + '\t' + str(int(x[2])+int(len(x[0]))) + '\t' + str(foundPromoters[item].begin) + '\t' 
            + str(foundPromoters[item].end) + '\t' + str(x[4]) + '\t' + "promoter"+ '\t' +str(x[5]))  #file write of promoters results
            if str(x[5]) == '0\n':
	            d_promoters_0MM[str(x[0])] +=1
            elif str(x[5]) == '1\n':
	            d_promoters_1MM[str(x[0])] +=1
            elif str(x[5]) == '2\n':
                d_promoters_2MM[str(x[0])] +=1
            elif str(x[5]) == '3\n':
                d_promoters_3MM[str(x[0])] +=1
            elif str(x[5]) == '4\n':
                d_promoters_4MM[str(x[0])] +=1
            elif str(x[5]) == '5\n':
                d_promoters_5MM[str(x[0])] +=1
            elif str(x[5]) == '6\n':
                d_promoters_6MM[str(x[0])] +=1
            elif str(x[5]) == '7\n':
                d_promoters_7MM[str(x[0])] +=1
            elif str(x[5]) == '8\n':
                d_promoters_8MM[str(x[0])] +=1
            elif str(x[5]) == '9\n':
                d_promoters_9MM[str(x[0])] +=1
            elif str(x[5]) == '10\n':
                d_promoters_10MM[str(x[0])] +=1
    
    for item in range(0,len(foundCtcf)):
        if foundCtcf[item].data == str(x[1]):
            outFile.write(str(x[3]) + '\t' + str(foundCtcf[item].data) + '\t' + str(x[2]) + '\t' + str(int(x[2])+int(len(x[0]))) + '\t' + str(foundCtcf[item].begin) + '\t' 
            + str(foundCtcf[item].end) + '\t' + str(x[4]) + '\t' + "ctcf"+ '\t' +str(x[5]))  #file write of ctcf results
            if str(x[5]) == '0\n':
	            d_ctcf_0MM[str(x[0])] +=1
            elif str(x[5]) == '1\n':
	            d_ctcf_1MM[str(x[0])] +=1
            elif str(x[5]) == '2\n':
                d_ctcf_2MM[str(x[0])] +=1
            elif str(x[5]) == '3\n':
                d_ctcf_3MM[str(x[0])] +=1
            elif str(x[5]) == '4\n':
                d_ctcf_4MM[str(x[0])] +=1
            elif str(x[5]) == '5\n':
                d_ctcf_5MM[str(x[0])] +=1
            elif str(x[5]) == '6\n':
                d_ctcf_6MM[str(x[0])] +=1
            elif str(x[5]) == '7\n':
                d_ctcf_7MM[str(x[0])] +=1
            elif str(x[5]) == '8\n':
                d_ctcf_8MM[str(x[0])] +=1
            elif str(x[5]) == '9\n':
                d_ctcf_9MM[str(x[0])] +=1
            elif str(x[5]) == '10\n':
                d_ctcf_10MM[str(x[0])] +=1
    
    for item in range(0,len(foundDnase)):
        if foundDnase[item].data == str(x[1]):
            outFile.write(str(x[3]) + '\t' + str(foundDnase[item].data) + '\t' + str(x[2]) + '\t' + str(int(x[2])+int(len(x[0]))) + '\t' + str(foundDnase[item].begin) + '\t' 
            + str(foundDnase[item].end) + '\t' + str(x[4]) + '\t' + "dnase"+ '\t' +str(x[5]))  #file write of dnase results
            if str(x[5]) == '0\n':
	            d_dnase_0MM[str(x[0])] +=1
            elif str(x[5]) == '1\n':
	            d_dnase_1MM[str(x[0])] +=1
            elif str(x[5]) == '2\n':
                d_dnase_2MM[str(x[0])] +=1
            elif str(x[5]) == '3\n':
                d_dnase_3MM[str(x[0])] +=1
            elif str(x[5]) == '4\n':
                d_dnase_4MM[str(x[0])] +=1
            elif str(x[5]) == '5\n':
                d_dnase_5MM[str(x[0])] +=1
            elif str(x[5]) == '6\n':
                d_dnase_6MM[str(x[0])] +=1
            elif str(x[5]) == '7\n':
                d_dnase_7MM[str(x[0])] +=1
            elif str(x[5]) == '8\n':
                d_dnase_8MM[str(x[0])] +=1
            elif str(x[5]) == '9\n':
                d_dnase_9MM[str(x[0])] +=1
            elif str(x[5]) == '10\n':
                d_dnase_10MM[str(x[0])] +=1

print("ANNOTATION CYCLE TIME: %s seconds" % (time.time() - start_time))


inGuides.seek(0)

outfileExons = open(annotatedFile + '.ExonsCount.txt', 'w')
outfileIntrons = open(annotatedFile + '.IntronsCount.txt',"w")
outfilePromoters = open(annotatedFile + '.PromotersCount.txt',"w")
outfileCTCF = open(annotatedFile + '.CTCFCount.txt',"w")
outfileDNAse = open(annotatedFile + '.DNAseCount.txt',"w")


for line in inGuides:
    x = line.split('\n')

    outfileExons.write(str(x[0]) + '\t' + str(d_exons_0MM[str(x[0])]) + '\t' + str(d_exons_1MM[str(x[0])]) + '\t' + str(d_exons_2MM[str(x[0])]) + '\t' 
    + str(d_exons_3MM[str(x[0])]) + '\t' + str(d_exons_4MM[str(x[0])]) + '\t' +str(d_exons_5MM[str(x[0])]) + '\t' + str(d_exons_6MM[str(x[0])]) + '\t'
    + str(d_exons_7MM[str(x[0])]) + '\t' + str(d_exons_8MM[str(x[0])]) + '\t' +str(d_exons_9MM[str(x[0])]) + '\t' + str(d_exons_10MM[str(x[0])]) + '\n')

    outfileIntrons.write(str(x[0]) + '\t' + str(d_introns_0MM[str(x[0])]) + '\t' + str(d_introns_1MM[str(x[0])]) + '\t' + str(d_introns_2MM[str(x[0])]) + '\t' 
    + str(d_introns_3MM[str(x[0])]) + '\t' + str(d_introns_4MM[str(x[0])]) + '\t' +str(d_introns_5MM[str(x[0])]) + '\t' + str(d_introns_6MM[str(x[0])]) + '\t'
    + str(d_introns_7MM[str(x[0])]) + '\t' + str(d_introns_8MM[str(x[0])]) + '\t' +str(d_introns_9MM[str(x[0])]) + '\t' + str(d_introns_10MM[str(x[0])]) + '\n')

    outfilePromoters.write(str(x[0]) + '\t' + str(d_promoters_0MM[str(x[0])]) + '\t' + str(d_promoters_1MM[str(x[0])]) + '\t' + str(d_promoters_2MM[str(x[0])]) + '\t' 
    + str(d_promoters_3MM[str(x[0])]) + '\t' + str(d_promoters_4MM[str(x[0])]) + '\t' +str(d_promoters_5MM[str(x[0])]) + '\t' + str(d_promoters_6MM[str(x[0])]) + '\t'
    + str(d_promoters_7MM[str(x[0])]) + '\t' + str(d_promoters_8MM[str(x[0])]) + '\t' +str(d_promoters_9MM[str(x[0])]) + '\t' + str(d_promoters_10MM[str(x[0])]) + '\n')

    outfileCTCF.write(str(x[0]) + '\t' + str(d_ctcf_0MM[str(x[0])]) + '\t' + str(d_ctcf_1MM[str(x[0])]) + '\t' + str(d_ctcf_2MM[str(x[0])]) + '\t' 
    + str(d_ctcf_3MM[str(x[0])]) + '\t' + str(d_ctcf_4MM[str(x[0])]) + '\t' +str(d_ctcf_5MM[str(x[0])]) + '\t' + str(d_ctcf_6MM[str(x[0])]) + '\t'
    + str(d_ctcf_7MM[str(x[0])]) + '\t' + str(d_ctcf_8MM[str(x[0])]) + '\t' +str(d_ctcf_9MM[str(x[0])]) + '\t' + str(d_ctcf_10MM[str(x[0])]) + '\n')

    outfileDNAse.write(str(x[0]) + '\t' + str(d_dnase_0MM[str(x[0])]) + '\t' + str(d_dnase_1MM[str(x[0])]) + '\t' + str(d_dnase_2MM[str(x[0])]) + '\t' 
    + str(d_dnase_3MM[str(x[0])]) + '\t' + str(d_dnase_4MM[str(x[0])]) + '\t' +str(d_dnase_5MM[str(x[0])]) + '\t' + str(d_dnase_6MM[str(x[0])]) + '\t'
    + str(d_dnase_7MM[str(x[0])]) + '\t' + str(d_dnase_8MM[str(x[0])]) + '\t' +str(d_dnase_9MM[str(x[0])]) + '\t' + str(d_dnase_10MM[str(x[0])]) + '\n')

#!/usr/bin/env python

from intervaltree import Interval, IntervalTree
import sys
import time
import concurrent.futures


print("READING")

# exonsFile = sys.argv[1]                 #exonsfile path
# intronsFile = sys.argv[2]               #intronsfile path
guidesFile = sys.argv[6]
resultFile = sys.argv[7]  # resultsfile path
annotatedFile = sys.argv[8]  # annotatedfile path
# inExons = open(exonsFile, "r")          #exonsfile open
# inIntrons = open(intronsFile, "r")      #intronsfile open
inResult = open(resultFile, "r")  # resultfile open
outFile = open(annotatedFile + '.targets.txt', 'w')  # outfile open
inGuides = open(guidesFile, "r")


start_time = time.time()


print("CREATING TREES")

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

total_0MM = 0
total_1MM = 0
total_2MM = 0
total_3MM = 0
total_4MM = 0
total_5MM = 0
total_6MM = 0
total_7MM = 0
total_8MM = 0
total_9MM = 0
total_10MM = 0

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

if(sys.argv[1] != "no"):
    exonsFile = sys.argv[1]
    inExons = open(exonsFile, "r")
    for line in inExons:
        x = line.split('\t')
        exonsTree[int(x[1]):int(x[2])] = str(x[0])

if(sys.argv[2] != "no"):
    intronsFile = sys.argv[2]
    inIntrons = open(intronsFile, "r")
    for line in inIntrons:
        x = line.split('\t')
        intronsTree[int(x[1]):int(x[2])] = str(x[0])

if(sys.argv[3] != "no"):
    promotersFile = sys.argv[3]
    inPromoters = open(promotersFile, "r")
    for line in inPromoters:
        x = line.split('\t')
        promotersTree[int(x[1]):int(x[2])] = str(x[0])

if(sys.argv[4] != "no"):
    ctcfFile = sys.argv[4]
    inCtcf = open(ctcfFile, "r")
    for line in inCtcf:
        x = line.split('\t')
        ctcfTree[int(x[1]):int(x[2])] = str(x[0])

if(sys.argv[5] != "no"):
    dnaseFile = sys.argv[5]
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

print("ANNOTATING")


outFile.write(
    "#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge_Size\tAnnotation_Type\n")

next(inResult)
for line in inResult:
    x = line.split('\t')
    x[1] = str(x[1]).replace("-","")
    # exons interval query
    foundExons = sorted(exonsTree[int(x[4]):(int(x[4])+int(len(x[1]))+1)])
    # introns interval query
    foundIntrons = sorted(intronsTree[int(x[4]):(int(x[4])+int(len(x[1]))+1)])
    foundPromoters = sorted(
        promotersTree[int(x[4]):(int(x[4])+int(len(x[1]))+1)])
    foundCtcf = sorted(ctcfTree[int(x[4]):(int(x[4])+int(len(x[1]))+1)])
    foundDnase = sorted(dnaseTree[int(x[4]):(int(x[4])+int(len(x[1]))+1)])

    # check if line in result falls in a tree for a genetic notation and updates the relative count and the summary count
    for item in range(0, len(foundExons)):
        if foundExons[item].data == str(x[3]):
            # outFile.write(str(x[3]) + '\t' + str(foundExons[item].data) + '\t' + str(x[2]) + '\t' + str(int(x[2])+int(len(x[0]))) + '\t' + str(foundExons[item].begin) + '\t'+ str(foundExons[item].end) + '\t' + str(x[4]) + '\t' + "exon" + '\t' + str(x[6]))  # file write of exons results
            outFile.write(line.rstrip() + '\t' + "exon\n")
            if int(x[6]) == 0:
                d_exons_0MM[str(x[1])] += 1
            elif int(x[6]) == 1:
                d_exons_1MM[str(x[1])] += 1
            elif int(x[6]) == 2:
                d_exons_2MM[str(x[1])] += 1
            elif int(x[6]) == 3:
                d_exons_3MM[str(x[1])] += 1
            elif int(x[6]) == 4:
                d_exons_4MM[str(x[1])] += 1
            elif int(x[6]) == 5:
                d_exons_5MM[str(x[1])] += 1
            elif int(x[6]) == 6:
                d_exons_6MM[str(x[1])] += 1
            elif int(x[6]) == 7:
                d_exons_7MM[str(x[1])] += 1
            elif int(x[6]) == 8:
                d_exons_8MM[str(x[1])] += 1
            elif int(x[6]) == 9:
                d_exons_9MM[str(x[1])] += 1
            elif int(x[6]) == 10:
                d_exons_10MM[str(x[1])] += 1

    for item in range(0, len(foundIntrons)):
        if foundIntrons[item].data == str(x[3]):
            # outFile.write(str(x[3]) + '\t' + str(foundIntrons[item].data) + '\t' + str(x[2]) + '\t' + str(int(x[2])+int(len(x[1]))) + '\t' + str(foundIntrons[item].begin) + '\t'+ str(foundIntrons[item].end) + '\t' + str(x[4]) + '\t' + "intron" + '\t' + str(x[6]))  # file write of introns results
            outFile.write(line.rstrip() + '\t' + "intron\n")
            if int(x[6]) == 0:
                d_introns_0MM[str(x[1])] += 1
            elif int(x[6]) == 1:
                d_introns_1MM[str(x[1])] += 1
            elif int(x[6]) == 2:
                d_introns_2MM[str(x[1])] += 1
            elif int(x[6]) == 3:
                d_introns_3MM[str(x[1])] += 1
            elif int(x[6]) == 4:
                d_introns_4MM[str(x[1])] += 1
            elif int(x[6]) == 5:
                d_introns_5MM[str(x[1])] += 1
            elif int(x[6]) == 6:
                d_introns_6MM[str(x[1])] += 1
            elif int(x[6]) == 7:
                d_introns_7MM[str(x[1])] += 1
            elif int(x[6]) == 8:
                d_introns_8MM[str(x[1])] += 1
            elif int(x[6]) == 9:
                d_introns_9MM[str(x[1])] += 1
            elif int(x[6]) == 10:
                d_introns_10MM[str(x[1])] += 1

    for item in range(0, len(foundPromoters)):
        if foundPromoters[item].data == str(x[3]):
            # outFile.write(str(x[3]) + '\t' + str(foundPromoters[item].data) + '\t' + str(x[2]) + '\t' + str(int(x[2])+int(len(x[1]))) + '\t' + str(foundPromoters[item].begin) + '\t'+ str(foundPromoters[item].end) + '\t' + str(x[4]) + '\t' + "promoter" + '\t' + str(x[6]))  # file write of promoters results
            outFile.write(line.rstrip() + '\t' + "promoter\n")
            if int(x[6]) == 0:
                d_promoters_0MM[str(x[1])] += 1
            elif int(x[6]) == 1:
                d_promoters_1MM[str(x[1])] += 1
            elif int(x[6]) == 2:
                d_promoters_2MM[str(x[1])] += 1
            elif int(x[6]) == 3:
                d_promoters_3MM[str(x[1])] += 1
            elif int(x[6]) == 4:
                d_promoters_4MM[str(x[1])] += 1
            elif int(x[6]) == 5:
                d_promoters_5MM[str(x[1])] += 1
            elif int(x[6]) == 6:
                d_promoters_6MM[str(x[1])] += 1
            elif int(x[6]) == 7:
                d_promoters_7MM[str(x[1])] += 1
            elif int(x[6]) == 8:
                d_promoters_8MM[str(x[1])] += 1
            elif int(x[6]) == 9:
                d_promoters_9MM[str(x[1])] += 1
            elif int(x[6]) == 10:
                d_promoters_10MM[str(x[1])] += 1

    for item in range(0, len(foundCtcf)):
        if foundCtcf[item].data == str(x[3]):
            # outFile.write(str(x[3]) + '\t' + str(foundCtcf[item].data) + '\t' + str(x[2]) + '\t' + str(int(x[2])+int(len(x[1]))) + '\t' + str(foundCtcf[item].begin) + '\t'+ str(foundCtcf[item].end) + '\t' + str(x[4]) + '\t' + "ctcf" + '\t' + str(x[6]))  # file write of ctcf results
            outFile.write(line.rstrip() + '\t' + "ctcf\n")
            if int(x[6]) == 0:
                d_ctcf_0MM[str(x[1])] += 1
            elif int(x[6]) == 1:
                d_ctcf_1MM[str(x[1])] += 1
            elif int(x[6]) == 2:
                d_ctcf_2MM[str(x[1])] += 1
            elif int(x[6]) == 3:
                d_ctcf_3MM[str(x[1])] += 1
            elif int(x[6]) == 4:
                d_ctcf_4MM[str(x[1])] += 1
            elif int(x[6]) == 5:
                d_ctcf_5MM[str(x[1])] += 1
            elif int(x[6]) == 6:
                d_ctcf_6MM[str(x[1])] += 1
            elif int(x[6]) == 7:
                d_ctcf_7MM[str(x[1])] += 1
            elif int(x[6]) == 8:
                d_ctcf_8MM[str(x[1])] += 1
            elif int(x[6]) == 9:
                d_ctcf_9MM[str(x[1])] += 1
            elif int(x[6]) == 10:
                d_ctcf_10MM[str(x[1])] += 1

    for item in range(0, len(foundDnase)):
        if foundDnase[item].data == str(x[3]):
            # outFile.write(str(x[3]) + '\t' + str(foundDnase[item].data) + '\t' + str(x[2]) + '\t' + str(int(x[2])+int(len(x[1]))) + '\t' + str(foundDnase[item].begin) + '\t'+ str(foundDnase[item].end) + '\t' + str(x[4]) + '\t' + "dnase" + '\t' + str(x[6]))  # file write of dnase results
            outFile.write(line.rstrip() + '\t' + "dnase\n")
            if int(x[6]) == 0:
                d_dnase_0MM[str(x[1])] += 1
            elif int(x[6]) == 1:
                d_dnase_1MM[str(x[1])] += 1
            elif int(x[6]) == 2:
                d_dnase_2MM[str(x[1])] += 1
            elif int(x[6]) == 3:
                d_dnase_3MM[str(x[1])] += 1
            elif int(x[6]) == 4:
                d_dnase_4MM[str(x[1])] += 1
            elif int(x[6]) == 5:
                d_dnase_5MM[str(x[1])] += 1
            elif int(x[6]) == 6:
                d_dnase_6MM[str(x[1])] += 1
            elif int(x[6]) == 7:
                d_dnase_7MM[str(x[1])] += 1
            elif int(x[6]) == 8:
                d_dnase_8MM[str(x[1])] += 1
            elif int(x[6]) == 9:
                d_dnase_9MM[str(x[1])] += 1
            elif int(x[6]) == 10:
                d_dnase_10MM[str(x[1])] += 1

    if int(x[6]) == 0:
        total_0MM += 1
    elif int(x[6]) == 1:
        total_1MM += 1
    elif int(x[6]) == 2:
        total_2MM += 1
    elif int(x[6]) == 3:
        total_3MM += 1
    elif int(x[6]) == 4:
        total_4MM += 1
    elif int(x[6]) == 5:
        total_5MM += 1
    elif int(x[6]) == 6:
        total_6MM += 1
    elif int(x[6]) == 7:
        total_7MM += 1
    elif int(x[6]) == 8:
        total_8MM += 1
    elif int(x[6]) == 9:
        total_9MM += 1
    elif int(x[6]) == 10:
        total_10MM += 1


print("ANNOTATION CYCLE TIME: %s seconds" % (time.time() - start_time))


inGuides.seek(0)

outfileExons = open(annotatedFile + '.ExonsCount.txt', "w")
outfileIntrons = open(annotatedFile + '.IntronsCount.txt', "w")
outfilePromoters = open(annotatedFile + '.PromotersCount.txt', "w")
outfileCTCF = open(annotatedFile + '.CTCFCount.txt', "w")
outfileDNAse = open(annotatedFile + '.DNAseCount.txt', "w")
outfileTotal = open(annotatedFile + '.SummaryCount.txt', "w")

outfileTotal.write("Total_targets" + '\t' + str(total_0MM) + '\t' + str(total_1MM) + '\t' + str(total_2MM) + '\t'
                   + str(total_3MM) + '\t' + str(total_4MM) + '\t'
                   + str(total_5MM) + '\t' + str(total_6MM) + '\t'
                   + str(total_7MM) + '\t' + str(total_8MM) + '\t'
                   + str(total_9MM) + '\t' + str(total_10MM) + '\n')

outfileTotal.write("Total_exons" + '\t' + str(sum(d_exons_0MM.values())) + '\t' + str(sum(d_exons_1MM.values())) + '\t' + str(sum(d_exons_2MM.values())) + '\t'
                   + str(sum(d_exons_3MM.values())) + '\t' +
                   str(sum(d_exons_4MM.values())) + '\t'
                   + str(sum(d_exons_5MM.values())) + '\t' +
                   str(sum(d_exons_6MM.values())) + '\t'
                   + str(sum(d_exons_7MM.values())) + '\t' +
                   str(sum(d_exons_8MM.values())) + '\t'
                   + str(sum(d_exons_9MM.values())) + '\t' + str(sum(d_exons_10MM.values())) + '\n')

outfileTotal.write("Total_introns" + '\t' + str(sum(d_introns_0MM.values())) + '\t' + str(sum(d_introns_1MM.values())) + '\t' + str(sum(d_introns_2MM.values())) + '\t'
                   + str(sum(d_introns_3MM.values())) + '\t' +
                   str(sum(d_introns_4MM.values())) + '\t'
                   + str(sum(d_introns_5MM.values())) + '\t' +
                   str(sum(d_introns_6MM.values())) + '\t'
                   + str(sum(d_introns_7MM.values())) + '\t' +
                   str(sum(d_introns_8MM.values())) + '\t'
                   + str(sum(d_introns_9MM.values())) + '\t' + str(sum(d_introns_10MM.values())) + '\n')

outfileTotal.write("Total_ctcf" + '\t' + str(sum(d_ctcf_0MM.values())) + '\t' + str(sum(d_ctcf_1MM.values())) + '\t' + str(sum(d_ctcf_2MM.values())) + '\t'
                   + str(sum(d_ctcf_3MM.values())) + '\t' +
                   str(sum(d_ctcf_4MM.values())) + '\t'
                   + str(sum(d_ctcf_5MM.values())) + '\t' +
                   str(sum(d_ctcf_6MM.values())) + '\t'
                   + str(sum(d_ctcf_7MM.values())) + '\t' +
                   str(sum(d_ctcf_8MM.values())) + '\t'
                   + str(sum(d_ctcf_9MM.values())) + '\t' + str(sum(d_ctcf_10MM.values())) + '\n')

outfileTotal.write("Total_dnase" + '\t' + str(sum(d_dnase_0MM.values())) + '\t' + str(sum(d_dnase_1MM.values())) + '\t' + str(sum(d_dnase_2MM.values())) + '\t'
                   + str(sum(d_dnase_3MM.values())) + '\t' +
                   str(sum(d_dnase_4MM.values())) + '\t'
                   + str(sum(d_dnase_5MM.values())) + '\t' +
                   str(sum(d_dnase_6MM.values())) + '\t'
                   + str(sum(d_dnase_7MM.values())) + '\t' +
                   str(sum(d_dnase_8MM.values())) + '\t'
                   + str(sum(d_dnase_9MM.values())) + '\t' + str(sum(d_dnase_10MM.values())) + '\n')

outfileTotal.write("Total_promoters" + '\t' + str(sum(d_promoters_0MM.values())) + '\t' + str(sum(d_promoters_1MM.values())) + '\t' + str(sum(d_promoters_2MM.values())) + '\t'
                   + str(sum(d_promoters_3MM.values())) + '\t' +
                   str(sum(d_promoters_4MM.values())) + '\t'
                   + str(sum(d_promoters_5MM.values())) + '\t' +
                   str(sum(d_promoters_6MM.values())) + '\t'
                   + str(sum(d_promoters_7MM.values())) + '\t' +
                   str(sum(d_promoters_8MM.values())) + '\t'
                   + str(sum(d_promoters_9MM.values())) + '\t' + str(sum(d_promoters_10MM.values())) + '\n')

for line in inGuides:
    x = line.split('\n')

    outfileExons.write(str(x[0]) + '\t' + str(d_exons_0MM[str(x[0])]) + '\t' + str(d_exons_1MM[str(x[0])]) + '\t' + str(d_exons_2MM[str(x[0])]) + '\t'
                       + str(d_exons_3MM[str(x[0])]) + '\t' +
                       str(d_exons_4MM[str(x[0])]) + '\t'
                       + str(d_exons_5MM[str(x[0])]) + '\t' +
                       str(d_exons_6MM[str(x[0])]) + '\t'
                       + str(d_exons_7MM[str(x[0])]) + '\t' + str(d_exons_8MM[str(x[0])]) + '\t' + str(d_exons_9MM[str(x[0])]) + '\t' + str(d_exons_10MM[str(x[0])]) + '\n')

    outfileIntrons.write(str(x[0]) + '\t' + str(d_introns_0MM[str(x[0])]) + '\t' + str(d_introns_1MM[str(x[0])]) + '\t' + str(d_introns_2MM[str(x[0])]) + '\t'
                         + str(d_introns_3MM[str(x[0])]) + '\t' + str(d_introns_4MM[str(x[0])]) + '\t' +
                         str(d_introns_5MM[str(x[0])]) + '\t' +
                         str(d_introns_6MM[str(x[0])]) + '\t'
                         + str(d_introns_7MM[str(x[0])]) + '\t' + str(d_introns_8MM[str(x[0])]) + '\t' + str(d_introns_9MM[str(x[0])]) + '\t' + str(d_introns_10MM[str(x[0])]) + '\n')

    outfilePromoters.write(str(x[0]) + '\t' + str(d_promoters_0MM[str(x[0])]) + '\t' + str(d_promoters_1MM[str(x[0])]) + '\t' + str(d_promoters_2MM[str(x[0])]) + '\t'
                           + str(d_promoters_3MM[str(x[0])]) + '\t' +
                           str(d_promoters_4MM[str(x[0])]) + '\t'
                           + str(d_promoters_5MM[str(x[0])]) + '\t' +
                           str(d_promoters_6MM[str(x[0])]) + '\t'
                           + str(d_promoters_7MM[str(x[0])]) + '\t' + str(d_promoters_8MM[str(x[0])]) + '\t' + str(d_promoters_9MM[str(x[0])]) + '\t' + str(d_promoters_10MM[str(x[0])]) + '\n')

    outfileCTCF.write(str(x[0]) + '\t' + str(d_ctcf_0MM[str(x[0])]) + '\t' + str(d_ctcf_1MM[str(x[0])]) + '\t' + str(d_ctcf_2MM[str(x[0])]) + '\t'
                      + str(d_ctcf_3MM[str(x[0])]) + '\t' +
                      str(d_ctcf_4MM[str(x[0])]) + '\t'
                      + str(d_ctcf_5MM[str(x[0])]) + '\t' +
                      str(d_ctcf_6MM[str(x[0])]) + '\t'
                      + str(d_ctcf_7MM[str(x[0])]) + '\t' + str(d_ctcf_8MM[str(x[0])]) + '\t' + str(d_ctcf_9MM[str(x[0])]) + '\t' + str(d_ctcf_10MM[str(x[0])]) + '\n')

    outfileDNAse.write(str(x[0]) + '\t' + str(d_dnase_0MM[str(x[0])]) + '\t' + str(d_dnase_1MM[str(x[0])]) + '\t' + str(d_dnase_2MM[str(x[0])]) + '\t'
                       + str(d_dnase_3MM[str(x[0])]) + '\t' +
                       str(d_dnase_4MM[str(x[0])]) + '\t'
                       + str(d_dnase_5MM[str(x[0])]) + '\t' +
                       str(d_dnase_6MM[str(x[0])]) + '\t'
                       + str(d_dnase_7MM[str(x[0])]) + '\t' + str(d_dnase_8MM[str(x[0])]) + '\t' + str(d_dnase_9MM[str(x[0])]) + '\t' + str(d_dnase_10MM[str(x[0])]) + '\n')

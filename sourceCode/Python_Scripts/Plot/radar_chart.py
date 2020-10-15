#!/usr/bin/env python

# Libraries
from operator import itemgetter
from os.path import isfile, join
from os import listdir
import os
import warnings
import glob
from itertools import islice
import sys
import numpy as np
import scipy.spatial.distance as sp
from math import pi
import pandas as pd
from matplotlib import patches as mpatches
from matplotlib import pyplot as plt
import math
import matplotlib
import random
# matplotlib.use("TkAgg")
matplotlib.use('Agg')

warnings.filterwarnings("ignore")

# argv 1 is Guide
# argv 2 is mm value
# argv 3 is annotation file (new version)
# argv 4 is extended profile (can be 'no' for web server)
# argv 5 is second summary for barplot (optional, 'no' if not given in input)
# argv 6 is gecko Summary file
# argv 7 is for web server: create png instead of pdf (-ws) (optional)
# argv 8 is for sample/pop/superpop report (-sample HG001/EUR/TSI) with this option, the annotation file is .sample_annotation.GUIDE.sample.txt
# TODO if argv8 is selected, from argv6 get the directory containing the summary of the specific sample for gecko -> line 618
# NOTE al momento lo script funziona con gecko.annotation.summary; la comparison con i sample la fa su annotation.summary. Quando l'analisi per sample sarà fatta, la comparison
# con i sampla dovrà essere fatta sul file specifico
plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["figure.figsize"] = [100, 80]
SIZE_GECKO = 123411  # NOTE modify if new gecko annotations are done
SIZE_GECKO = 111671
random.seed(a=None, version=2)
guide = sys.argv[1]
mm = int(sys.argv[2])
# annotation_var = sys.argv[3] #sostituito da lettura automatica del file summary

# utilizzato per calcolare la 'bontà' della guida usata in comparazione con il suo set (opzionale quello delle popolazioni)
# annotation delle guide ref, tutte le guide usate nella ricerca
annotation_ref = sys.argv[3]
inAnnotation_ref = open(annotation_ref, 'r').readlines()
# utilizzato per motif logo
extended_profile = sys.argv[4]
# annotation delle popolazioni, un file per guida
annotation_population = sys.argv[5]
inAnnotation_population = open(annotation_population, 'r').readlines()


# gecko_summary = sys.argv[6]
# gecko_exon = sys.argv[7]
# gecko_intron = sys.argv[8]
# gecko_promoter = sys.argv[9]
# gecko_dnase = sys.argv[10]
# gecko_ctcf = sys.argv[11]

web_server = False
population = False
file_extension = 'pdf'

if '-ws' in sys.argv[:]:
    web_server = True
if web_server:
    file_extension = 'png'
if '-pop' is sys.argv[:]:
    population = True

dictRef = {}
dictPop = {}
extractGuide = ''
pop = ''
populationSet = set()
annotationSet = set()
populationDataFrame = pd.DataFrame()

for line in inAnnotation_ref:
    split = line.lower().strip().split('\t')
    if 'Summary' in line:
        extractGuide = line.strip().split('_')[1]
        dictRef[extractGuide+'_REF'] = {}
    else:
        # dictRef[extractGuide+'_REF'][split[0]] = int(split[mm+1])
        dictRef[extractGuide+'_REF'][split[0]] = random.randint(200, 123456)
        annotationSet.add(split[0].lower())
        annotationSet.add(split[0].lower()+'_rank')

dictRef.pop('Total_REF', None)
# print(dictRef)

referenceDataFrame = pd.DataFrame.from_dict(dictRef, orient='index')
# print(referenceDataFrame)


for col in referenceDataFrame:
    referenceDataFrame[str(
        col)+'_rank'] = referenceDataFrame[col].rank(method='first', pct=True)

for line in inAnnotation_population:
    split = line.lower().strip().split('\t')
    if 'Guide' in line:
        extractGuide = line.strip().split('_')[2]
        dictPop[extractGuide] = {}
    elif 'Summary' in line and 'Guide' not in line:
        pop = line.strip().split('_')[1]
        populationSet.add(pop.upper())
        dictPop[extractGuide][pop] = {}
    else:
        # dictPop[extractGuide][pop][split[0]] = int(split[mm+1])
        dictPop[extractGuide][pop][split[0]] = random.randint(200, 123456)

for key in dictPop:
    dictPop[key].pop('Total', None)

populationSet.remove('TOTAL')

tempDict = {}
for elem in populationSet:
    for key in dictPop:
        tempDict[key+'_'+elem] = {}
        tempDict[key+'_'+elem] = dictPop[key][elem]
    tempDataFrame = pd.DataFrame.from_dict(tempDict, orient='index')
    for col in tempDataFrame:
        tempDataFrame[str(
            col)+'_rank'] = tempDataFrame[col].rank(method='first', pct=True)
    populationDataFrame = pd.concat([populationDataFrame, tempDataFrame])
    tempDict = {}

# print(referenceDataFrame)
# print(populationDataFrame)

# print(annotationSet)

totalDataFrame = pd.concat([referenceDataFrame, populationDataFrame])
# totalDataFrame.to_csv('total.tsv', sep='\t')

totalDataFrame = totalDataFrame.T

populationSet = sorted(populationSet)
populationSet.append('REF')
annotationSet = sorted(annotationSet)
# print(totalDataFrame)
# Create data for radarchart
dataForChart = {'Populations': populationSet}
dataForTable = []

for annot in annotationSet:
    if 'rank' in annot:
        dataForChart[annot] = []

# print(dataForChart)
# data_for_table_df = [dictRef['targets']]
# rows = ['General']
# listForChart = []
templist = []
for pop in populationSet:
    for col in totalDataFrame:
        if pop in col and guide in col:
            for annot in annotationSet:
                templist.append(totalDataFrame.loc[annot, col])
                if 'rank' in annot:
                    dataForChart[annot].append(totalDataFrame.loc[annot, col])
            dataForTable.append(templist)
            templist = []

# print(dataForChart)
dataFrameForChart = pd.DataFrame.from_dict(dataForChart, orient='index')
dataFrameForChart = dataFrameForChart.T
# print(dataFrameForChart)
# number of variable
categories = list(dataFrameForChart)[1:]
N = len(categories)

# # We are going to plot the first line of the data frame.
# # But we need to repeat the first value to close the circular graph:
# values = dataFrameForChart.loc[0].drop('').values.flatten().tolist()
# values += values[:1]

# # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]

# # Initialise the spider plot
# ax = plt.subplot(111, polar=True)
ax = plt.subplot(1, 1, 1, polar=True)
# plt.title('RADAR CHART')

labels = list(dataFrameForChart.columns.values[1:])
# Draw one axe per variable + add labels labels yet
plt.xticks(angles[:-1], labels, color='black', size=18)

for label, rot in zip(ax.get_xticklabels(), angles):
    if (rot == 0):
        label.set_horizontalalignment("center")
    if (rot > 0):
        label.set_horizontalalignment("left")
    if (rot > 3):
        label.set_horizontalalignment("center")
    if (rot > 4):
        label.set_horizontalalignment("right")

# # offset posizione y-axis
ax.set_theta_offset(pi / 2)
ax.set_theta_direction(-1)

# # Draw ylabels
ax.set_rlabel_position(0)
plt.yticks([0, 0.25, 0.50, 0.75], ["0", "0.25",
                                   "0.50", "0.75"], color="black", size=12)
plt.ylim(0, 1)

# listPop = list(populationSet)
for count, pop in enumerate(populationSet):
    values = dataFrameForChart.loc[count].drop(
        'Populations').values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, linewidth=0.1,
            linestyle='solid', label=str(pop))
    ax.fill(angles, values, alpha=0.1)

# # Ind2
# values = dataFrameForChart.loc[1].drop('Populations').values.flatten().tolist()
# values += values[:1]
# ax.plot(angles, values, linewidth=1, linestyle='solid', label=str(listPop[1]))
# ax.fill(angles, values, 'r', alpha=0.1)

# Add legend
# plt.legend(bbox_to_anchor=(0.1, 0.1))


# # Plot data
# ax.plot(angles, values, linewidth=1, linestyle='solid')

# # Fill area
# ax.fill(angles, values, 'b', alpha=0.1)

# columns = ('Position', '# Targets')

# print(dataForTable)
# print(dataForChart)

numpy_array = np.array(dataForTable)
transpose = numpy_array.T
transpose_list = transpose.tolist()
# print(transpose_list)

tableAnnotation = []
for annot in annotationSet:
    if 'rank' not in annot:
        tableAnnotation.append(annot)

# print(tableAnnotation)

# Create table plot
plt.subplot(1, 1, 1)
table = plt.table(cellText=transpose_list, rowLabels=annotationSet,
                  colLabels=annotationSet, loc='center', colWidths=[0.30 for x in annotationSet])
table.auto_set_font_size(False)
table.set_fontsize(18)
table.scale(1, 3)
plt.axis('off')


plt.savefig("summary_single_guide_" + str(guide) + "_"+str(mm) +
            "mm" + "." + file_extension, format=file_extension)

# # Calculate motif logo
# if extended_profile != 'no':
#     # lists for data storing and analysis
#     guidesExtendedProfile = []
#     # reading extendend profile to obtain results over mismatches counts
#     with open(extended_profile, 'r') as inGuidesProfileExtended:
#         for line in inGuidesProfileExtended:
#             if ">" + guide in line:
#                 # print(line)
#                 next(inGuidesProfileExtended)
#                 # line=inGuidesProfileExtended.readline()
#                 for ciao in range(0, uppermm+1):
#                     line = inGuidesProfileExtended.readline().strip()
#                     count = 0
#                     x = line.split('\t')
#                     # read X Mismatch line, exclude last column
#                     guidesExtendedProfile.append(tuple(x[1:-1]))
#                     # guidesExtendedProfile.append((x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9],
#                     #                             x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20]))
#                     for line in inGuidesProfileExtended:
#                         if count < 6:
#                             line = line.rstrip()
#                             x = line.split('\t')
#                             # y = str(x[20]).split('\n')
#                             guidesExtendedProfile.append(tuple(x[1:]))
#                             # guidesExtendedProfile.append((x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9],
#                             #                             x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20]))
#                             count += 1
#                         else:
#                             break
#                 break
#     arrayguidesExtendedProfile = np.array(guidesExtendedProfile, dtype=int)
#     shape_size = len(guidesExtendedProfile[0])
#     arrayguidesExtendedProfile.shape = (7*((uppermm-0)+1), shape_size)

#     # Plot motif logo
#     datacount = arrayguidesExtendedProfile[mm*7] / \
#         (max(arrayguidesExtendedProfile[mm*7]))
#     data = np.array(datacount, dtype=float)
#     data = np.around(data, decimals=1)
#     data.shape = (1, len(datacount))

#     # string = guide[0:20]
#     string = guide.replace('N', '')
#     strArray = np.array([list(string)])

#     A = arrayguidesExtendedProfile[mm*7+1] / \
#         (max(arrayguidesExtendedProfile[mm*7]))
#     C = arrayguidesExtendedProfile[mm*7+2] / \
#         (max(arrayguidesExtendedProfile[mm*7]))
#     G = arrayguidesExtendedProfile[mm*7+3] / \
#         (max(arrayguidesExtendedProfile[mm*7]))
#     T = arrayguidesExtendedProfile[mm*7+4] / \
#         (max(arrayguidesExtendedProfile[mm*7]))
#     DNA = arrayguidesExtendedProfile[mm*7+5] / \
#         (max(arrayguidesExtendedProfile[mm*7]))
#     RNA = arrayguidesExtendedProfile[mm*7+6] / \
#         (max(arrayguidesExtendedProfile[mm*7]))

#     # the x locations for the groups
#     ind = np.arange(0, len(string), 1) + 0.15
#     width = 0.7  # the width of the bars: can also be len(x) sequence

#     motif = plt.subplot(2, 1, 2, frameon=False)
#     # motif=plt.subplot(1,1,1)
#     p1 = plt.bar(ind, A, width, color='red', align='edge')
#     p2 = plt.bar(ind, C, width, color='blue', bottom=A, align='edge')
#     p3 = plt.bar(ind, G, width, color='green', bottom=A+C, align='edge')
#     p4 = plt.bar(ind, T, width, color='purple', bottom=C+G+A, align='edge')
#     p5 = plt.bar(ind, DNA, width, color='magenta',
#                  bottom=C+G+A+T, align='edge')
#     p6 = plt.bar(ind, RNA, width, color='gold',
#                  bottom=C+G+A+T+DNA, align='edge')
#     plt.xlim(0, len(string))
#     plt.xticks([])

#     plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]),
#                ('A', 'C', 'G', 'T', 'D', 'R'), fontsize=18)

#     table = plt.table(cellText=strArray, loc='bottom',
#                       cellLoc='center', rowLoc='bottom')
#     table.auto_set_font_size(False)
#     table.set_fontsize(18)
#     table.scale(1, 1.6)
#     table.xticks = ([])
#     table.yticks = ([])

# plt.suptitle(str(mm)+" Mismatches",
#              horizontalalignment='center', color='black', size=25)

# plt.tight_layout()
# plt.subplots_adjust(top=0.90, bottom=0.07, left=0.05,
#                     right=0.99, wspace=0.12)

# plt.savefig("summary_single_guide_" + str(guide) +
#             "_"+str(mm) + "mm" + "." + file_extension, format=file_extension)

# elif radarchart_sample:  # Radarchart for sample
#     if gecko_summary != 'no':  # NOTE sta usando il summary gecko totale, non quello diviso per sample
#         summary_two = []
#         with open(annotation_var, 'r') as ann_var:
#             for line in ann_var:
#                 if sample_name in line:
#                     found = True
#                     continue
#                 if found:
#                     if 'Summary_' in line:
#                         break
#                     summary_two.append(line.strip())

#         # Save annotation shared between gecko and user annotations
#         common_annotations = dict()
#         for i in summary_two:
#             if i.split('\t')[0].lower() in ['targets', 'ctcf', 'dnasei', 'exon', 'intron', 'promoter']:
#                 # {'targets': [3, 7, 25, 474, 6872, 0, 0, 0, 0, 0], 'ctcf_binding_site': [0, 0, 3, 61, 659, 0, 0, 0, 0, 0],...}
#                 common_annotations[i.split('\t')[0].lower()] = [
#                     int(x) for x in i.split('\t')[1:]]

#         ################ NEW GECKO COMPARISON #########################
#         with open(gecko_summary, 'r') as ann_var:
#             # Create the complete matrix (row = guide, col = annotation, cell = number of targets with that guide in that annotation with mm mismatches)
#             current_guide = ''
#             for line in ann_var:
#                 line = line.strip().split('\t')
#                 if line[0] != current_guide and '_' in line[0]:
#                     complete_matrix.append([line[0].split('_')[1]])
#                     current_guide = line[0]
#                     continue
#                 # Add to complete matrix only common annotations
#                 if line[0].lower() in common_annotations:
#                     complete_matrix[-1].append(line[mm + 1])
#         # NOTE -1 to eliminate the 'Summary_ Total' row that is always the last one when the matrix will be sorted
#         num_guides = len(complete_matrix) - 1
#         min_diff = sys.maxsize
#         pos_in_gecko = SIZE_GECKO
#         cluster_pos = SIZE_GECKO

#         # Order the 1, 2 , 3 ... col and get input guide position
#         # ['targets','CTCF', 'DNAseI', 'Exons', 'Introns', 'Promoters']    #NOTE lexicograph order (same as summary_two)
#         annotations = sorted_common_annotations
#         min_diff = sys.maxsize
#         pos_in_gecko = SIZE_GECKO
#         cluster_pos = SIZE_GECKO

#         for pos, annotation in enumerate(annotations):
#             # int(summary_two[pos].split('\t')[mm + 1])
#             target_value = common_annotations[annotation][mm]
#             complete_matrix.sort(key=lambda x: int(
#                 x[pos + 1]))  # +1 pechè 0 è nome guida
#             for pos_best, row in enumerate(complete_matrix):

#                 current_diff = abs(int(row[pos+1]) - target_value)
#                 if current_diff > min_diff:
#                     pos_in_gecko = cluster_pos  # pos - 1
#                     break
#                 elif current_diff < min_diff:
#                     cluster_pos = pos_best

#                 min_diff = current_diff

#             table_dict[annotation] = [
#                 round(pos_in_gecko / SIZE_GECKO, 2), target_value]
#             min_diff = sys.maxsize
#             pos_in_gecko = SIZE_GECKO
#             cluster_pos = SIZE_GECKO

#         # OLD GECKO COMPARISON WITH SUMMARY DIVIDED IN FILES ################# #NOTE posso riusare questo codice per calc sample
#             # con gecko, basta cambiare line[0] in 'qualsiasi cosa'
#             # e mm+1 con mm, quando aggiungo a_l

#         gecko_summary = gecko_summary[:gecko_summary.rfind(
#             '/')] + '/Summary_Gecko/' + sample_name
#         # annotations = ['CTCF', 'DNAseI', 'Exons', 'Introns', 'Promoters']    #NOTE lexicograph order (same as summary_two)
#         annotations = ['CTCF', 'dnasei', 'exon', 'intron', 'promoter']
#         # for ann_pos, annotation in enumerate(annotations):
#         #     target_value = int(summary_two[ann_pos + 1].split('\t')[mm + 1])     #first +1 due to targets elelment , second +1 due to name of annotation in splitted list
#         #     #Load file
#         #     with open(gecko_summary + '/gecko.' + annotation + '.Count.txt', 'r') as ge: <- no
#         # with open (gecko_summary + '/Summary_Gecko/' + sample_name + '/' + annotation + '.count_gecko.txt') as ge: <- this
#         #         a_l = []
#         #         for line in ge:
#         #             line = line.strip().split('\t')
#         #             a_l.append([line[0], int(line[mm + 1])]) <- no       #+1 due to name of annotation in line.split
#         # a_l.append(['0', int(line[mm])]) <- this       # '0' is just placeholder
#         #     a_l.sort(key=lambda x: int(x[1]))
#         #     #Calculate min difference
#         #     min_diff = sys.maxsize
#         #     pos_in_gecko = SIZE_GECKO
#         #     cluster_pos = SIZE_GECKO
#         #     for pos, i in enumerate(a_l):
#         #         current_diff = abs(i[1] - target_value)
#         #         if current_diff > min_diff:
#         #             pos_in_gecko = cluster_pos #pos - 1
#         #             break
#         #         elif current_diff < min_diff:
#         #             cluster_pos = pos

#         #         min_diff = current_diff

#         #     #table_dict[annotation.lower()] = [round(pos_in_gecko / SIZE_GECKO, 2), target_value]
#         #     table_dict[annotation] = [round(pos_in_gecko / SIZE_GECKO, 2), target_value]

#         #     #Calculate position on Total targets value
#         #     a_l = []
#         #     with open(gecko_summary + '/gecko.reference.profile.xls') as gecko_pro:
#         #         header = gecko_pro.readline().strip().split('\t')
#         #         column_mm = header.index(str(mm) +  'MM')
#         #         for line in gecko_pro:
#         #             line = line.strip().split('\t')
#         #             a_l.append([line[0], int(line[column_mm])])
#         #     a_l.sort(key=lambda x: int(x[1]))
#         #     #Calculate min difference
#         #     target_value = int(summary_two[0].split('\t')[mm + 1])
#         #     min_diff = sys.maxsize
#         #     pos_in_gecko = SIZE_GECKO
#         #     cluster_pos = SIZE_GECKO
#         #     for pos, i in enumerate(a_l):
#         #         current_diff = abs(i[1] - target_value)
#         #         if current_diff > min_diff:
#         #             pos_in_gecko = cluster_pos#pos - 1
#         #             break
#         #         elif current_diff < min_diff:
#         #             cluster_pos = pos
#         #         min_diff = current_diff

#         #     table_dict['targets'] = [round(pos_in_gecko / SIZE_GECKO,2), target_value]
#     else:

#         # summary_two contains counting of total, change variable name to summary_total
#         summary_total = summary_two.copy()
#         # Extract counting for total
#         summary_two = []
#         with open(annotation_var, 'r') as ann_var:
#             for line in ann_var:
#                 if sample_name in line:
#                     found = True
#                     continue
#                 if found:
#                     if 'Summary_' in line:
#                         break
#                     summary_two.append(line.strip())
#         # Calculate distances
#         # print('num target: targets', summary_total[0].split('\t')[mm + 1])
#         if int(summary_total[0].split('\t')[mm + 1]) > 0:
#             table_dict = {'targets': [round(1 - round(int(summary_two[0].split('\t')[mm + 1]) / int(
#                 summary_total[0].split('\t')[mm + 1]), 2), 2), summary_two[0].split('\t')[mm + 1]]}
#         else:
#             table_dict = {'targets': [0.0, summary_two[0].split('\t')[mm + 1]]}
#         for pos, annotation in enumerate(summary_two[1:], 1):
#             annotation = annotation.split('\t')
#             # print('numtarget:', annotation,summary_total[pos].split('\t')[mm + 1] )
#             if int(summary_total[pos].split('\t')[mm + 1]) > 0:
#                 table_dict[annotation[0]] = [round(1 - round(int(annotation[mm + 1]) / int(
#                     summary_total[pos].split('\t')[mm + 1]), 2), 2), annotation[mm + 1]]
#             else:
#                 table_dict[annotation[0]] = [0.0, annotation[mm + 1]]

#     # Create data for radarchart
#     data_for_df = {'group': ['A'], 'General': table_dict['targets'][0]}
#     data_for_table_df = [table_dict['targets']]
#     rows = ['General']
#     for elem in table_dict:
#         if elem == 'targets':
#             continue
#         data_for_df[elem] = table_dict[elem][0]
#         data_for_table_df.append(table_dict[elem])
#         rows.append(elem)

#     dataFrameForChart = pd.DataFrame(data_for_df)

#     # number of variable
#     categories = list(dataFrameForChart)[1:]
#     N = len(categories)

#     # We are going to plot the first line of the data frame.
#     # But we need to repeat the first value to close the circular graph:
#     values = dataFrameForChart.loc[0].drop('group').values.flatten().tolist()
#     values += values[:1]

#     # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
#     angles = [n / float(N) * 2 * pi for n in range(N)]
#     angles += angles[:1]

#     # Initialise the spider plot
#     ax = plt.subplot(1, 2, 1, polar=True)
#     # ax=plt.subplot(1, 1, 1, polar=True)
#     # plt.title('RADAR CHART')

#     labels = list(dataFrameForChart.columns.values[1:])
#     # Draw one axe per variable + add labels labels yet
#     plt.xticks(angles[:-1], labels, color='black', size=18)
#     for label, rot in zip(ax.get_xticklabels(), angles):
#         if (rot == 0):
#             label.set_horizontalalignment("center")
#         if (rot > 0):
#             label.set_horizontalalignment("left")
#         if (rot > 3):
#             label.set_horizontalalignment("center")
#         if (rot > 4):
#             label.set_horizontalalignment("right")

#     # offset posizione y-axis
#     ax.set_theta_offset(pi / 2)
#     ax.set_theta_direction(-1)

#     # Draw ylabels
#     ax.set_rlabel_position(0)
#     plt.yticks([0, 0.25, 0.50, 0.75, 1], ["0", "0.25",
#                                           "0.50", "0.75"], color="black", size=18)
#     plt.ylim(0, 1)

#     # Plot data
#     ax.plot(angles, values, linewidth=1, linestyle='solid')

#     # Fill area
#     ax.fill(angles, values, 'b', alpha=0.1)

#     columns = ('Position', '# Targets')

#     # Create table plot
#     plt.subplot(1, 2, 2)
#     table = plt.table(cellText=data_for_table_df, rowLabels=rows,
#                       colLabels=columns, loc='center', colWidths=[0.35 for x in columns])
#     table.auto_set_font_size(False)
#     table.set_fontsize(18)
#     table.scale(1, 3)
#     plt.axis('off')

#     plt.suptitle(str(mm)+" Mismatches" + ' - ' + sample_name,
#                  horizontalalignment='center', color='black', size=25)

#     plt.tight_layout()
#     plt.subplots_adjust(top=0.90, bottom=0.07, left=0.10,
#                         right=0.99, wspace=0.50)

#     plt.savefig("summary_single_guide_" + str(sample_name) +
#                 "_" + guide_name_for_sample + '_' + str(mm) + "mm" + "." + file_extension, format=file_extension)

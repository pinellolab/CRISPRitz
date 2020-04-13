#!/usr/bin/env python

# Libraries
import math
import matplotlib
matplotlib.use("TkAgg")
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import pandas as pd
from math import pi
import scipy.spatial.distance as sp
import numpy as np
import sys
from itertools import islice
import glob
import warnings
import os
from os import listdir
from os.path import isfile, join
warnings.filterwarnings("ignore")

from operator import itemgetter

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
SIZE_GECKO = 123411 #NOTE modify if new gecko annotations are done
SIZE_GECKO = 111671
guide = sys.argv[1]
mm = int(sys.argv[2])
annotation_var = sys.argv[3]
extended_profile = sys.argv[4]

annotation_ref = sys.argv[5]
gecko_summary = sys.argv[6]
# gecko_exon = sys.argv[7]
# gecko_intron = sys.argv[8]
# gecko_promoter = sys.argv[9]
# gecko_dnase = sys.argv[10]
# gecko_ctcf = sys.argv[11]

if '-ws' in sys.argv[:]:
    web_server = True
else:
    web_server = False
file_extension ='pdf'
if web_server:
    file_extension = 'png'

radarchart_sample = False
guide_name_for_sample = ''
if '-sample' in sys.argv[:]:
    
    sample_name = sys.argv.index('-sample')
    sample_name = sys.argv[sample_name + 1]
    radarchart_sample = True
    guide_name_for_sample = guide
    guide = 'no'


summary_two = []        #NOTE summary_two is the summary of the variant genome
summary_one = []        #NOTE summary one is the summary of the reference genome
found = False

if guide == 'no':
    search_term = 'Total'
else:
    search_term = guide
with open(annotation_var, 'r') as ann_var:
    for line in ann_var:
        if search_term in line:
            found = True
            continue
        if found:
            if 'Summary_' in line:
                break
            summary_two.append(line.strip())

found = False          
if annotation_ref != 'no':
    with open(annotation_ref, 'r') as ann_ref:
    
        for line in ann_ref:
            if search_term in line:
                found = True
                continue
            if found:
                if 'Summary_' in line:
                    break
                summary_one.append(line.strip())

uppermm = mm
if summary_one  and summary_two :           #DO BARPLOT
    inSummaryCountOne = summary_one #open(summary_one, "r")
    inSummaryCountTwo = summary_two #open(summary_two, "r")
    arraySummaryCountOne = []
    arraySummaryCountTwo = []
    
    x_ticks_labels = []
    n_annotation = 0
    mms_total = 0      #len of the array exon   0   0   0   0   0   0   0, taken from [1:]
    for line in inSummaryCountOne:
        x = line.strip().split('\t')
        arraySummaryCountOne.append(tuple(x[1:]))
        x_ticks_labels.append(x[0])
        n_annotation = n_annotation + 1
        mms_total = len(x[1:])
    

    for line in inSummaryCountTwo:
        x = line.strip().split('\t')
        arraySummaryCountTwo.append(tuple(x[1:]))
    arraySummaryCountOne = np.array(arraySummaryCountOne, dtype=int)
    arraySummaryCountOne.shape = (n_annotation, mms_total)

    arraySummaryCountTwo = np.array(arraySummaryCountTwo, dtype=int)
    arraySummaryCountTwo.shape = (n_annotation, mms_total)

    percentageGain = []
    for row in range(0, n_annotation):
        for col in range(0, uppermm+1):
            res = max((arraySummaryCountOne[row, col]/arraySummaryCountTwo[row, col]),
                      (arraySummaryCountTwo[row, col]/arraySummaryCountOne[row, col]))
            percentageGain.append(res)
    percentageGain = np.array(percentageGain, dtype=float)
    percentageGain.shape = (n_annotation, uppermm+1)

    intergenicGainOne = (np.sum(arraySummaryCountOne, axis=0))
    intergenicGainTwo = (np.sum(arraySummaryCountTwo, axis=0))
    intergenicGainOne.shape = (1, mms_total)
    intergenicGainTwo.shape = (1, mms_total)
    intergenicGainOne = intergenicGainOne - arraySummaryCountOne[0, :]
    intergenicGainTwo = intergenicGainTwo - arraySummaryCountTwo[0, :]

    # p1 = plt.bar(1, percentageGain[0, uppermm-1])
    # p2 = plt.bar(2, percentageGain[1, uppermm-1])
    # p3 = plt.bar(3, percentageGain[2, uppermm-1])
    # p4 = plt.bar(4, percentageGain[3, uppermm-1])
    # p5 = plt.bar(5, percentageGain[4, uppermm-1])
    # p6 = plt.bar(6, percentageGain[5, uppermm-1])
    # p7 = plt.bar(7, max((intergenicGainOne[0, uppermm-1]/intergenicGainTwo[0, uppermm-1]),
    #                     (intergenicGainTwo[0, uppermm-1]/intergenicGainOne[0, uppermm-1])))

    ind = np.arange(0, n_annotation, 1)
    # print('math ceil', np.arange(0, max(arraySummaryCountTwo[:, uppermm]) + math.ceil(max(
    #     arraySummaryCountTwo[:, uppermm])/10), 
    #     math.ceil(max(arraySummaryCountTwo[:, uppermm])/5)))  #BUG if max(arraySummaryCountTwo[:, uppermm], then arange (0,0,0)
    no_result = False
    try:
        y_range = np.arange(0, max(arraySummaryCountTwo[:, uppermm]) + math.ceil(max(
            arraySummaryCountTwo[:, uppermm])/10), math.ceil(max(arraySummaryCountTwo[:, uppermm])/5))
    except:
        y_range = np.arange(0,1,1)
        no_result = True
    width = 0.5

    p1 = plt.bar(
        ind, arraySummaryCountOne[:, uppermm], width, color='#67a9cf', align='edge')
    p2 = plt.bar(ind, (arraySummaryCountTwo[:, uppermm]-arraySummaryCountOne[:, uppermm]),
                 width, bottom=arraySummaryCountOne[:, uppermm], color='#ef8a62', align='edge')

    plt.legend((p1[0], p2[0]), ('Reference Genome',
                                'Enriched Genome'), fontsize=30)

    # plt.xlim(0, len(string))
    # plt.set_ylim([1, 1.5])
    plt.title('Relative Increase Enriched/Reference Genome with ' +
              str(uppermm) + ' Mismatches', size=25)
    # plt.xlabel('Annotations')
    if no_result:
        plt.annotate('No targets found with ' + str(mm)  + ' mismatches', [1.35,0], size = 22) #NOTE with 0-mm print only the mm pdf; 1.35 modificare se cambia la str
    else:
        for k in range(0, n_annotation):
            plt.annotate('%.2fx' % percentageGain[k, uppermm], [
                        k+0.05, arraySummaryCountTwo[k, uppermm]+(max(arraySummaryCountTwo[:, uppermm])/100)], size=22)
        # plt.ylim([0, max(arraySummaryCountTwo[:, uppermm])+2000],size=25)
    mx_len_ann = 0
    for i in x_ticks_labels:
        if len(i) > mx_len_ann:
            mx_len_ann = len(i)
    if mx_len_ann > 15:
        size_x = 10
    else:
        size_x = 22
    plt.xticks(ind+0.25, x_ticks_labels, size=size_x)
    plt.yticks(y_range, size=22)

    plt.tight_layout()
    plt.subplots_adjust(top=0.95, bottom=0.06, left=0.1, right=0.99)
    if guide != 'no':
        plt.savefig("summary_histogram_" + str(guide) + '_' + str(uppermm) + 
                    "mm" + "." + file_extension, format=file_extension)
    else:
        plt.savefig("summary_histogram_" + str(uppermm) + 
                    "mm" + "." + file_extension, format=file_extension)
    # if guide != 'no':
    #     plt.savefig("summary_histogram_" + str(guide) + '_' + str(uppermm) + 
    #                 "mm" + ".png", format="png")
    # else:
    #     plt.savefig("summary_histogram_" + str(uppermm) + 
    #                 "mm" + ".png", format="png")
    


check_annotation_name = []
common_annotations = dict()     #Save annotation shared between gecko and user annotations
sorted_common_annotations = []
for i in summary_two:
    check_annotation_name.append(i.split('\t')[0].lower())
    if i.split('\t')[0].lower() in ['targets', 'ctcf', 'dnasei', 'exon', 'intron', 'promoter']:
        common_annotations[i.split('\t')[0].lower()] = [int(x) for x in i.split('\t')[1:]] ##{'targets': [3, 7, 25, 474, 6872, 0, 0, 0, 0, 0], 'ctcf_binding_site': [0, 0, 3, 61, 659, 0, 0, 0, 0, 0],...}
        sorted_common_annotations.append(i.split('\t')[0].lower())  # NOTE targets always first
if set(check_annotation_name) != set(['targets', 'ctcf', 'dnasei', 'exon', 'intron', 'promoter']) and gecko_summary != "no":     #NOTE if gecko annotations are updated, update the list adding the new annotations 
    print('Warning! Option \'-gecko\' can be used on files annotated with \'ctcf\', \'dnasei\', \'exon\', \'intron\', \'promoter\'')
    print('Using -gecko only on common annotations')
    # gecko_summary = 'no'

complete_matrix = []  
table_dict = dict()    
if guide != 'no':       #DO RADAR CHART
    if gecko_summary != 'no':       
        ########################### With New Gecko Summary #######################################
        with open(gecko_summary, 'r') as ann_var:
            #Create the complete matrix (row = guide, col = annotation, cell = number of targets with that guide in that annotation with mm mismatches)
            current_guide = ''
            for line in ann_var:
                line = line.strip().split('\t')
                if line[0] != current_guide and '_' in line[0]:
                    complete_matrix.append([line[0].split('_')[1]])
                    current_guide = line[0]
                    continue
                if line[0].lower() in common_annotations:       #Add to complete matrix only common annotations
                    complete_matrix[-1].append(line[mm + 1])
    
        num_guides = len(complete_matrix) - 1   #NOTE -1 to eliminate the 'Summary_ Total' row that is always the last one when the matrix will be sorted
        min_diff = sys.maxsize
        pos_in_gecko = SIZE_GECKO
        cluster_pos = SIZE_GECKO
        #Order the 1, 2 , 3 ... col and get input guide position
        annotations = sorted_common_annotations #['targets','CTCF', 'DNAseI', 'Exons', 'Introns', 'Promoters']    #NOTE lexicograph order (same as summary_two)
        min_diff = sys.maxsize
        pos_in_gecko = SIZE_GECKO
        cluster_pos = SIZE_GECKO
        
        for pos, annotation in enumerate(annotations):
            target_value = common_annotations[annotation][mm] #int(summary_two[pos].split('\t')[mm + 1])
            complete_matrix.sort(key=lambda x: int(x[pos + 1 ]))     #+1 pechè 0 è nome guida
            for pos_best, row in enumerate(complete_matrix):
                current_diff = abs(int(row[pos+1]) - target_value)
                if current_diff > min_diff:
                    pos_in_gecko = cluster_pos #pos - 1
                    break
                elif current_diff < min_diff:
                    cluster_pos = pos_best
                min_diff = current_diff
            table_dict[annotation] = [round(pos_in_gecko / SIZE_GECKO, 2), target_value]
            min_diff = sys.maxsize
            pos_in_gecko = SIZE_GECKO
            cluster_pos = SIZE_GECKO
        
        ############## Old Gecko Counting, with SummaryCount.Exon etc files ###############

        # #Calculate position on Total targets value
        # a_l = []
        # with open(gecko_summary + '/gecko.reference.profile.xls') as gecko_pro:
        #     header = gecko_pro.readline().strip().split('\t')
        #     column_mm = header.index(str(mm) +  'MM')
        #     for line in gecko_pro:
        #         line = line.strip().split('\t')
        #         a_l.append([line[0], int(line[column_mm])])
        # a_l.sort(key=lambda x: int(x[1]))
        # #Calculate min difference
        # target_value = int(summary_two[0].split('\t')[mm + 1])
        # min_diff = sys.maxsize
        # pos_in_gecko = SIZE_GECKO
        # cluster_pos = SIZE_GECKO
        # for pos, i in enumerate(a_l):
        #     current_diff = abs(i[1] - target_value)
        #     if current_diff > min_diff:
        #         pos_in_gecko = cluster_pos#pos - 1
        #         break
        #     elif current_diff < min_diff:
        #         cluster_pos = pos
        #     min_diff = current_diff

        # table_dict['targets'] = [round(pos_in_gecko / SIZE_GECKO,2), target_value]
        # min_diff = sys.maxsize
        # pos_in_gecko = SIZE_GECKO
        # cluster_pos = SIZE_GECKO
        # for pos, i in enumerate(a_l):
        #     target_value = int(summary_two[ann_pos + 1].split('\t')[mm + 1]) 
        #     current_diff = abs(i[1] - target_value)
        #     if current_diff > min_diff:
        #         pos_in_gecko = cluster_pos #pos - 1
        #         break
        #     elif current_diff < min_diff:
        #         cluster_pos = pos

        #     min_diff = current_diff
        
        # table_dict[annotation.lower()] = [round(pos_in_gecko / SIZE_GECKO, 2), target_value]
        
        
        
        
        # ###################################################################
        # #do process for gecko -> annotations file before 08/12/2019
        # annotations = ['CTCF', 'DNAseI', 'Exons', 'Introns', 'Promoters']    #NOTE lexicograph order (same as summary_two)
        # for ann_pos, annotation in enumerate(annotations):
        #     target_value = int(summary_two[ann_pos + 1].split('\t')[mm + 1])     #first +1 due to targets elelment , second +1 due to name of annotation in splitted list
        #     #Load file
        #     with open(gecko_summary + '/gecko.' + annotation + '.Count.txt', 'r') as ge:
        #         a_l = []
        #         for line in ge:
        #             line = line.strip().split('\t')
        #             a_l.append([line[0], int(line[mm + 1])])        #+1 due to name of annotation in line.split
            
        #     a_l.sort(key=lambda x: int(x[1]))
            
        #     #Calculate min difference
        #     min_diff = sys.maxsize
        #     pos_in_gecko = SIZE_GECKO
        #     cluster_pos = SIZE_GECKO
        #     for pos, i in enumerate(a_l):
        #         current_diff = abs(i[1] - target_value)
        #         if current_diff > min_diff:
        #             pos_in_gecko = cluster_pos #pos - 1
        #             break
        #         elif current_diff < min_diff:
        #             cluster_pos = pos

        #         min_diff = current_diff
            
        #     # table_dict[annotation.lower()] = [round(pos_in_gecko / SIZE_GECKO, 2), target_value]
        #     table_dict[annotation] = [round(pos_in_gecko / SIZE_GECKO, 2), target_value]
        
        # #Calculate position on Total targets value
        # a_l = []
        # with open(gecko_summary + '/gecko.reference.profile.xls') as gecko_pro:
        #     header = gecko_pro.readline().strip().split('\t')
        #     column_mm = header.index(str(mm) +  'MM')
        #     for line in gecko_pro:
        #         line = line.strip().split('\t')
        #         a_l.append([line[0], int(line[column_mm])])
        # a_l.sort(key=lambda x: int(x[1]))
        # #Calculate min difference
        # target_value = int(summary_two[0].split('\t')[mm + 1])
        # min_diff = sys.maxsize
        # pos_in_gecko = SIZE_GECKO
        # cluster_pos = SIZE_GECKO
        # for pos, i in enumerate(a_l):
        #     current_diff = abs(i[1] - target_value)
        #     if current_diff > min_diff:
        #         pos_in_gecko = cluster_pos#pos - 1
        #         break
        #     elif current_diff < min_diff:
        #         cluster_pos = pos
        #     min_diff = current_diff

        # table_dict['targets'] = [round(pos_in_gecko / SIZE_GECKO,2), target_value]
            
            
    else:
        with open(annotation_var, 'r') as ann_var:
            #Create the complete matrix (row = guide, col = annotation, cell = number of targets with that guide in that annotation with mm mismatches)
            current_guide = ''
            for line in ann_var:
                line = line.strip().split('\t')
                if line[0] != current_guide and '-Summary_' in line[0]:
                    complete_matrix.append([line[0].split('_')[1]])
                    current_guide = line[0]
                    continue
                complete_matrix[-1].append(line[mm + 1])
        num_guides = len(complete_matrix) - 1   #NOTE -1 to eliminate the 'Summary_Total' row that is always the last one

        #Order the 1, 2 , 3 ... col and get input guide position
        for pos, annotation in enumerate(check_annotation_name):
            complete_matrix.sort(key=lambda x: int(x[pos +1 ]))
            for pos_best, row in enumerate(complete_matrix):
                if guide == row[0]:
                    table_dict[annotation] = [round(pos_best / num_guides, 2), row[pos + 1]]
                    break 

    #Create data for radarchart
    data_for_df = {'group': ['A'], 'General':table_dict['targets'][0]}
    data_for_table_df = [table_dict['targets']]
    rows = ['General']
    for elem in table_dict:
        if elem == 'targets':
            continue
        data_for_df[elem] = table_dict[elem][0]
        data_for_table_df.append(table_dict[elem])
        rows.append(elem)
    
    df = pd.DataFrame(data_for_df)

    # number of variable
    categories = list(df)[1:]
    N = len(categories)

    # We are going to plot the first line of the data frame.
    # But we need to repeat the first value to close the circular graph:
    values = df.loc[0].drop('group').values.flatten().tolist()
    values += values[:1]

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    # Initialise the spider plot
    ax = plt.subplot(2, 2, 1, polar=True)
    # ax=plt.subplot(1, 1, 1, polar=True)
    # plt.title('RADAR CHART')

    labels = list(df.columns.values[1:])
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

    # offset posizione y-axis
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    # Draw ylabels
    ax.set_rlabel_position(0)
    plt.yticks([0, 0.25, 0.50, 0.75, 1], ["0", "0.25",
                                        "0.50", "0.75"], color="black", size=18)
    plt.ylim(0, 1)

    # Plot data
    ax.plot(angles, values, linewidth=1, linestyle='solid')

    # Fill area
    ax.fill(angles, values, 'b', alpha=0.1)

    columns = ('Position', '# Targets')
    
    #Create table plot
    plt.subplot(2, 2, 2)
    table = plt.table(cellText=data_for_table_df, rowLabels=rows,colLabels=columns, loc='center', colWidths=[0.35 for x in columns])
    table.auto_set_font_size(False)
    table.set_fontsize(18)
    table.scale(1, 3)
    plt.axis('off')

    #Calculate motif logo
    if extended_profile != 'no':
        # lists for data storing and analysis
        guidesExtendedProfile = []
        # reading extendend profile to obtain results over mismatches counts
        with open(extended_profile, 'r') as inGuidesProfileExtended:
            for line in inGuidesProfileExtended:
                if ">" + guide in line:
                    # print(line)
                    next(inGuidesProfileExtended)
                    # line=inGuidesProfileExtended.readline()
                    for ciao in range(0, uppermm+1):
                        line = inGuidesProfileExtended.readline().strip()
                        count = 0
                        x = line.split('\t')
                        guidesExtendedProfile.append(tuple(x[1:-1]))    #read X Mismatch line, exclude last column
                        # guidesExtendedProfile.append((x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9],
                        #                             x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20]))
                        for line in inGuidesProfileExtended:
                            if count < 6:
                                line = line.rstrip()
                                x = line.split('\t')
                                #y = str(x[20]).split('\n')
                                guidesExtendedProfile.append(tuple(x[1:]))
                                # guidesExtendedProfile.append((x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9],
                                #                             x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20]))
                                count += 1
                            else:
                                break
                    break
        arrayguidesExtendedProfile = np.array(guidesExtendedProfile, dtype=int)
        shape_size = len(guidesExtendedProfile[0])
        arrayguidesExtendedProfile.shape = (7*((uppermm-0)+1), shape_size)

        #Plot motif logo
        datacount = arrayguidesExtendedProfile[mm*7] / \
                (max(arrayguidesExtendedProfile[mm*7]))
        data = np.array(datacount, dtype=float)
        data = np.around(data, decimals=1)
        data.shape = (1, len(datacount))

        # string = guide[0:20]
        string = guide.replace('N','')
        strArray = np.array([list(string)])

        A = arrayguidesExtendedProfile[mm*7+1] / \
            (max(arrayguidesExtendedProfile[mm*7]))
        C = arrayguidesExtendedProfile[mm*7+2] / \
            (max(arrayguidesExtendedProfile[mm*7]))
        G = arrayguidesExtendedProfile[mm*7+3] / \
            (max(arrayguidesExtendedProfile[mm*7]))
        T = arrayguidesExtendedProfile[mm*7+4] / \
            (max(arrayguidesExtendedProfile[mm*7]))
        DNA = arrayguidesExtendedProfile[mm*7+5] / \
            (max(arrayguidesExtendedProfile[mm*7]))
        RNA = arrayguidesExtendedProfile[mm*7+6] / \
            (max(arrayguidesExtendedProfile[mm*7]))

        # the x locations for the groups
        ind = np.arange(0, len(string), 1) + 0.15
        width = 0.7  # the width of the bars: can also be len(x) sequence

        motif = plt.subplot(2, 1, 2, frameon=False)
        # motif=plt.subplot(1,1,1)
        p1 = plt.bar(ind, A, width, color='red', align='edge')
        p2 = plt.bar(ind, C, width, color='blue', bottom=A, align='edge')
        p3 = plt.bar(ind, G, width, color='green', bottom=A+C, align='edge')
        p4 = plt.bar(ind, T, width, color='purple', bottom=C+G+A, align='edge')
        p5 = plt.bar(ind, DNA, width, color='magenta', bottom=C+G+A+T, align='edge')
        p6 = plt.bar(ind, RNA, width, color='gold', bottom=C+G+A+T+DNA, align='edge')
        plt.xlim(0, len(string))
        plt.xticks([])

        plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]),
                ('A', 'C', 'G', 'T', 'D', 'R'), fontsize=18)

        table = plt.table(cellText=strArray, loc='bottom',
                        cellLoc='center', rowLoc='bottom')
        table.auto_set_font_size(False)
        table.set_fontsize(18)
        table.scale(1, 1.6)
        table.xticks = ([])
        table.yticks = ([])
    
    plt.suptitle(str(mm)+" Mismatches",
                horizontalalignment='center', color='black', size=25)

    plt.tight_layout()
    plt.subplots_adjust(top=0.90, bottom=0.07, left=0.05,
                        right=0.99, wspace=0.12)

    plt.savefig("summary_single_guide_" + str(guide) +
                "_"+str(mm) + "mm" + "." + file_extension, format=file_extension)

elif radarchart_sample:   #Radarchart for sample
    if gecko_summary != 'no':       #NOTE sta usando il summary gecko totale, non quello diviso per sample
        summary_two = []
        with open(annotation_var, 'r') as ann_var:
            for line in ann_var:
                if sample_name in line:
                    found = True
                    continue
                if found:
                    if 'Summary_' in line:
                        break
                    summary_two.append(line.strip())
        
        common_annotations = dict()     #Save annotation shared between gecko and user annotations
        for i in summary_two:
            if i.split('\t')[0].lower() in ['targets', 'ctcf', 'dnasei', 'exon', 'intron', 'promoter']:
                common_annotations[i.split('\t')[0].lower()] = [int(x) for x in i.split('\t')[1:]] ##{'targets': [3, 7, 25, 474, 6872, 0, 0, 0, 0, 0], 'ctcf_binding_site': [0, 0, 3, 61, 659, 0, 0, 0, 0, 0],...}

        ################ NEW GECKO COMPARISON #########################
        with open(gecko_summary, 'r') as ann_var:
            #Create the complete matrix (row = guide, col = annotation, cell = number of targets with that guide in that annotation with mm mismatches)
            current_guide = ''
            for line in ann_var:
                line = line.strip().split('\t')
                if line[0] != current_guide and '_' in line[0]:
                    complete_matrix.append([line[0].split('_')[1]])
                    current_guide = line[0]
                    continue
                if line[0].lower() in common_annotations:       #Add to complete matrix only common annotations
                    complete_matrix[-1].append(line[mm + 1])
        num_guides = len(complete_matrix) - 1   #NOTE -1 to eliminate the 'Summary_ Total' row that is always the last one when the matrix will be sorted
        min_diff = sys.maxsize
        pos_in_gecko = SIZE_GECKO
        cluster_pos = SIZE_GECKO
        
        #Order the 1, 2 , 3 ... col and get input guide position
        annotations = sorted_common_annotations#['targets','CTCF', 'DNAseI', 'Exons', 'Introns', 'Promoters']    #NOTE lexicograph order (same as summary_two)
        min_diff = sys.maxsize
        pos_in_gecko = SIZE_GECKO
        cluster_pos = SIZE_GECKO

        for pos, annotation in enumerate(annotations):
            target_value = common_annotations[annotation][mm] #int(summary_two[pos].split('\t')[mm + 1])
            complete_matrix.sort(key=lambda x: int(x[pos + 1 ]))     #+1 pechè 0 è nome guida
            for pos_best, row in enumerate(complete_matrix):
                
                current_diff = abs(int(row[pos+1]) - target_value)
                if current_diff > min_diff:
                    pos_in_gecko = cluster_pos #pos - 1
                    break
                elif current_diff < min_diff:
                    cluster_pos = pos_best

                min_diff = current_diff
            
            table_dict[annotation] = [round(pos_in_gecko / SIZE_GECKO, 2), target_value]
            min_diff = sys.maxsize
            pos_in_gecko = SIZE_GECKO
            cluster_pos = SIZE_GECKO





        ########### OLD GECKO COMPARISON WITH SUMMARY DIVIDED IN FILES ################# #NOTE posso riusare questo codice per calc sample 
                                                                                        # con gecko, basta cambiare line[0] in 'qualsiasi cosa'
                                                                                        # e mm+1 con mm, quando aggiungo a_l
        
        gecko_summary = gecko_summary[:gecko_summary.rfind('/')] + '/Summary_Gecko/' + sample_name
        # annotations = ['CTCF', 'DNAseI', 'Exons', 'Introns', 'Promoters']    #NOTE lexicograph order (same as summary_two)
        annotations = ['CTCF', 'dnasei', 'exon', 'intron', 'promoter']
        # for ann_pos, annotation in enumerate(annotations):
        #     target_value = int(summary_two[ann_pos + 1].split('\t')[mm + 1])     #first +1 due to targets elelment , second +1 due to name of annotation in splitted list
        #     #Load file
        #     with open(gecko_summary + '/gecko.' + annotation + '.Count.txt', 'r') as ge: <- no
            #with open (gecko_summary + '/Summary_Gecko/' + sample_name + '/' + annotation + '.count_gecko.txt') as ge: <- this
        #         a_l = []
        #         for line in ge:
        #             line = line.strip().split('\t')
        #             a_l.append([line[0], int(line[mm + 1])]) <- no       #+1 due to name of annotation in line.split
                        #a_l.append(['0', int(line[mm])]) <- this       # '0' is just placeholder
        #     a_l.sort(key=lambda x: int(x[1]))
        #     #Calculate min difference
        #     min_diff = sys.maxsize
        #     pos_in_gecko = SIZE_GECKO
        #     cluster_pos = SIZE_GECKO
        #     for pos, i in enumerate(a_l):
        #         current_diff = abs(i[1] - target_value)
        #         if current_diff > min_diff:
        #             pos_in_gecko = cluster_pos #pos - 1
        #             break
        #         elif current_diff < min_diff:
        #             cluster_pos = pos

        #         min_diff = current_diff
            
        #     #table_dict[annotation.lower()] = [round(pos_in_gecko / SIZE_GECKO, 2), target_value]
        #     table_dict[annotation] = [round(pos_in_gecko / SIZE_GECKO, 2), target_value]
            
        #     #Calculate position on Total targets value
        #     a_l = []
        #     with open(gecko_summary + '/gecko.reference.profile.xls') as gecko_pro:
        #         header = gecko_pro.readline().strip().split('\t')
        #         column_mm = header.index(str(mm) +  'MM')
        #         for line in gecko_pro:
        #             line = line.strip().split('\t')
        #             a_l.append([line[0], int(line[column_mm])])
        #     a_l.sort(key=lambda x: int(x[1]))
        #     #Calculate min difference
        #     target_value = int(summary_two[0].split('\t')[mm + 1])
        #     min_diff = sys.maxsize
        #     pos_in_gecko = SIZE_GECKO
        #     cluster_pos = SIZE_GECKO
        #     for pos, i in enumerate(a_l):
        #         current_diff = abs(i[1] - target_value)
        #         if current_diff > min_diff:
        #             pos_in_gecko = cluster_pos#pos - 1
        #             break
        #         elif current_diff < min_diff:
        #             cluster_pos = pos
        #         min_diff = current_diff

        #     table_dict['targets'] = [round(pos_in_gecko / SIZE_GECKO,2), target_value]
    else:

        #summary_two contains counting of total, change variable name to summary_total
        summary_total = summary_two.copy()
        #Extract counting for total
        summary_two = []
        with open(annotation_var, 'r') as ann_var:
            for line in ann_var:
                if sample_name in line:
                    found = True
                    continue
                if found:
                    if 'Summary_' in line:
                        break
                    summary_two.append(line.strip())
        #Calculate distances
        # print('num target: targets', summary_total[0].split('\t')[mm + 1])
        if int(summary_total[0].split('\t')[mm + 1]) > 0:
            table_dict= {'targets':[round(1 - round(int(summary_two[0].split('\t')[mm + 1]) / int(summary_total[0].split('\t')[mm + 1]),2),2), summary_two[0].split('\t')[mm + 1]]}
        else:
            table_dict= {'targets': [0.0,summary_two[0].split('\t')[mm + 1] ]}
        for pos, annotation in enumerate(summary_two[1:], 1):
            annotation = annotation.split('\t')
            # print('numtarget:', annotation,summary_total[pos].split('\t')[mm + 1] )
            if int(summary_total[pos].split('\t')[mm + 1]) > 0:
                table_dict[annotation[0]]=[round(1 - round(int(annotation[mm + 1]) / int(summary_total[pos].split('\t')[mm + 1]),2),2), annotation[mm + 1]]
            else:
                table_dict[annotation[0]] = [0.0, annotation[mm + 1]]
    
    #Create data for radarchart
    data_for_df = {'group': ['A'], 'General':table_dict['targets'][0]}
    data_for_table_df = [table_dict['targets']]
    rows = ['General']
    for elem in table_dict:
        if elem == 'targets':
            continue
        data_for_df[elem] = table_dict[elem][0]
        data_for_table_df.append(table_dict[elem])
        rows.append(elem)
    
    df = pd.DataFrame(data_for_df)

    # number of variable
    categories = list(df)[1:]
    N = len(categories)

    # We are going to plot the first line of the data frame.
    # But we need to repeat the first value to close the circular graph:
    values = df.loc[0].drop('group').values.flatten().tolist()
    values += values[:1]

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    # Initialise the spider plot
    ax = plt.subplot(1, 2, 1, polar=True)
    # ax=plt.subplot(1, 1, 1, polar=True)
    # plt.title('RADAR CHART')

    labels = list(df.columns.values[1:])
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

    # offset posizione y-axis
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    # Draw ylabels
    ax.set_rlabel_position(0)
    plt.yticks([0, 0.25, 0.50, 0.75, 1], ["0", "0.25",
                                        "0.50", "0.75"], color="black", size=18)
    plt.ylim(0, 1)

    # Plot data
    ax.plot(angles, values, linewidth=1, linestyle='solid')

    # Fill area
    ax.fill(angles, values, 'b', alpha=0.1)

    columns = ('Position', '# Targets')
    
    #Create table plot
    plt.subplot(1, 2, 2)
    table = plt.table(cellText=data_for_table_df, rowLabels=rows,colLabels=columns, loc='center', colWidths=[0.35 for x in columns])
    table.auto_set_font_size(False)
    table.set_fontsize(18)
    table.scale(1, 3)
    plt.axis('off')
    


    plt.suptitle(str(mm)+" Mismatches" + ' - ' + sample_name,
                horizontalalignment='center', color='black', size=25)

    plt.tight_layout()
    plt.subplots_adjust(top=0.90, bottom=0.07, left=0.10,
                        right=0.99, wspace=0.50)

    plt.savefig("summary_single_guide_" + str(sample_name) +
                "_"+ guide_name_for_sample + '_' + str(mm) + "mm" + "." + file_extension, format=file_extension)
    


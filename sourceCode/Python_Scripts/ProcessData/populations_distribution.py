#!/usr/bin/env python


'''
Generate barplot with population distribution from the PopulationDistribution.txt file. Create a barplot for each total value, dividing
into 0 bulge, 1 bulge, 2 bulge etc
Example of PopulationDistribution file (with 1 bulge max):
-Summary_CCATCGGTGGCCGTTTGCCCNNN
EAS     0,0     0,0     0,0     0,0     0,4     0,22    0,0     0,0     0,0     0,0
EUR     0,0     0,0     0,0     0,0     2,1     0,28    0,0     0,0     0,0     0,0
AFR     0,0     1,0     0,1     0,0     1,6     0,42    0,0     0,0     0,0     0,0
AMR     0,0     0,0     0,0     0,0     2,4     0,30    0,0     0,0     0,0     0,0
SAS     0,0     0,0     0,0     0,0     1,1     0,34    0,0     0,0     0,0     0,0
-Summary_GAGTCCGAGCAGAAGAAGAANNN
EAS     0,0     0,0     0,0     1,1     1,18    0,147   0,0     0,0     0,0     0,0
EUR     0,0     0,0     0,0     1,1     3,20    0,128   0,0     0,0     0,0     0,0
AFR     0,0     0,0     0,0     1,1     3,31    0,253   0,0     0,0     0,0     0,0
AMR     0,0     0,0     0,0     1,1     3,26    0,146   0,0     0,0     0,0     0,0
SAS     0,0     0,0     0,0     2,1     3,21    0,153   0,0     0,0     0,0     0,0
'''

# argv 1 is jobid.PopulationDistribution.txt
# argv 2 is total value
# argv 3 is Guide
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
import colorsys
import matplotlib.colors as mc
from matplotlib.legend_handler import  HandlerTuple


def adjust_lightness(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

barplot_values = dict()         #barplot_values -> EAS -> [1,2]
total = int(sys.argv[2])
guide = sys.argv[3]
max_value = 0
number_bars = 0
previous_bar = []
result_name = (sys.argv)[1].split('.PopulationDistribution')[0]
with open(sys.argv[1]) as summary:
    for line in summary:
        if guide in line:
            while True:
                try:
                    line = next(summary).strip()
                except:
                    break
                if '-Summary_' in line:
                    break
                line = line.split('\t')
                number_bars = len(line[total + 1].split(','))
                barplot_values[line[0]] = [int(x) for x in line[total + 1].split(',')]       #line = EAS 0,7 1,2 5,3 10,11                
                value = sum([int(x) for x in line[total + 1].split(',')] )
                previous_bar.append(0)
                if value > max_value:
                    max_value = value

ind = np.arange(0, len(barplot_values.keys()), 1)       #[0 1 2 3 4 5] #NOTE 0 is REFERENCE
no_result = False
try:
    y_range = np.arange(0, max_value + math.ceil(max_value/10), math.ceil(max_value/5))
except:
    y_range = np.arange(0,1,1)
    no_result = True
width = 0.5

population_color = ['purple', 'orange', 'green', 'blue', 'red']
if 'REFERENCE' in barplot_values.keys():
    population_color.insert(0, 'grey')
all_bar = []
for i in range(number_bars):     #For 0 bulge, 1 bulge, 2 bulge ...
    current_bar = [x[i] for x in barplot_values.values()] 
    all_bar.append(plt.bar(ind, current_bar, width, color=[adjust_lightness(x, 1 + i*0.3) if current_bar[pos] != 0 else 'white' for pos, x in enumerate(population_color)], align='edge', bottom = previous_bar, edgecolor = 'black'))
    previous_bar = [ x + previous_bar[k] for k, x in enumerate(current_bar)]

# p1 = plt.bar(ind, barplot_values.values(), width, color=population_color, align='edge')       #color = '#67a9cf'

# p2 = plt.bar(ind, barplot_values.values() , width, color=[adjust_lightness(x, 1.3) for x in population_color], align='edge', bottom = list(barplot_values.values()))

legend_labels = []
handles_color = []
for x in range(min(number_bars, total+1)):
    if x == 0:
        legend_labels.append(str(total) + ' MM')            #got Total mms and 0 bulges
        handles_color.append((all_bar[0][:]))
    else:
        if (total - x) < 0:                                 #To avoid negative numbers on MM values
            legend_labels.append('0 MM + ' + str(x) + ' B')
            handles_color.append((all_bar[x][:]))
        else:
            legend_labels.append(str(total - x ) + ' MM + ' + str(x) + ' B')
            handles_color.append((all_bar[x][:]))
legend_labels.reverse()
# handles_color = [(x[:]) for x in all_bar]
handles_color.reverse()
plt.legend(handles_color, legend_labels, fontsize=15, handlelength = 10,handler_map={tuple: HandlerTuple(ndivide=None)})
# first param is for the colored rectangles of legens, second parameter for labels, handlelength is size of rectangles, handlermap is for grouping different colors in single label
#[(first bar color, second bar color, ...), (first bar light color, second bar light color,...)]
plt.title('Targets found in each Superpopulation - ' + str(total) + ' Mismatches(MM) + Bulges(B)', size=23)

if no_result:
    plt.annotate('No targets found with ' + str(total)  + ' mismatches + bulges', (2.5,0) ,size = 22, ha = 'center', va = 'center')  #2.5 is x position
    sys.exit()


plt.xticks(ind+0.25, barplot_values.keys(), size=25)

size_y_ticks = 22
digits = int(math.log10(max_value))+1
if digits > 5:  #Reduce dimension of y label because it can exceed plot size
    size_y_ticks = max(16, size_y_ticks - (2*(digits-5)))

plt.yticks(y_range, size=size_y_ticks)

plt.tight_layout()
plt.subplots_adjust(top=0.95, bottom=0.06, left=0.1, right=0.99)
plt.savefig(result_name + ".populations_distribution_" + guide + '_' + str(total) + "total" + ".png", format='png')
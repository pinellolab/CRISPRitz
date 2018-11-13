#!/usr/bin/python3

# Libraries
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from math import pi
import numpy as np
import sys
from statsmodels.distributions.empirical_distribution import ECDF
from itertools import islice
import glob
import warnings
warnings.filterwarnings("ignore")


if len(sys.argv)<6:
    print ("example of call: python radar_chart.py guidesProfileFile guidesExtendedProfileFile exonsCountFile intronsCountFile promotersCountFile DNAseCountFile CTCFCountFile guidetoanalyze missmatch")
    exit()


#lettura file                      
guidesProfileFile = sys.argv[1]
guidesExtendedProfileFile = sys.argv[2]
exonsCountFile = sys.argv[3]
intronsCountFile = sys.argv[4]
promotersCountFile = sys.argv[5]
DNAseCountFile = sys.argv[6]
CTCFCountFile = sys.argv[7]
guide = str(sys.argv[8])
if len(sys.argv[9])==1:
   missmatch=int(sys.argv[9])
   lowermm=0
   uppermm=missmatch
else:
   missmatch=sys.argv[9].split('-')
   lowermm=int(missmatch[0])
   uppermm=int(missmatch[1])
geckoProfile=sys.argv[10]
geckoExonsCount=sys.argv[11]
geckoIntronsCount=sys.argv[12]
geckoPromotersCount=sys.argv[13]
geckoDNAseCount=sys.argv[14]
geckoCTCFCount=sys.argv[15]
  

#apertura file 
inGuidesProfile = open(guidesProfileFile, "r")
inGuidesProfileExtended = open(guidesExtendedProfileFile, "r")   
inExonsCountFile = open(exonsCountFile, "r")
inIntronsCountFile = open(intronsCountFile, "r")  
inPromotersCountFile = open(promotersCountFile, "r")
inDNAseCountFile = open(DNAseCountFile, "r")
inCTCFCountFile = open(CTCFCountFile, "r")

#lists for data storing and analysis
guidesExtendedProfile = []

#global counting for annotation types
profileMissmatchGlobal = []
exonsMissmatchGlobal = []
intronsMissmatchGlobal = []
promotersMissmatchGlobal = []
DNAseMissmatchGlobal = []
CTCFMissmatchGlobal = []

#Empirical distribution function for annotation types
ecdfProfile = []
ecdfExons = []
ecdfIntrons = []
ecdfPromoters = []
ecdfDNAse = []
ecdfCTCF = []

#percentile calculated for annotation types
percentileProfile = []
percentileExons = []
percentileIntrons = []
percentilePromoters = []
percentileDNAse = []
percentileCTCF = []


#reading extendend profile to obtain results over mismatches counts
for line in inGuidesProfileExtended:
   if ">"+guide in line:
      #print(line)
      next(inGuidesProfileExtended)
      #line=inGuidesProfileExtended.readline()
      for ciao in range(0,uppermm+1):
         line=inGuidesProfileExtended.readline()
         #print(line)
         count = 0
         x=line.split('\t')
         guidesExtendedProfile.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20]))
         for line in inGuidesProfileExtended:
            if count<4:
               #print("count ",count,line)
               x = line.split('\t')
               y = str(x[20]).split('\n')
               guidesExtendedProfile.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],y[0]))
               count+=1
            else:
               break
      break


arrayguidesExtendedProfile = np.array(guidesExtendedProfile,dtype=int)
arrayguidesExtendedProfile.shape = (5*((uppermm-0)+1),20)


#reading profile file to obtain results for every mismatch count in the general profile
next(inGuidesProfile)
for line in inGuidesProfile:
    x = line.split('\t')
    profileMissmatchGlobal.append((x[26],x[27],x[28],x[29],x[30],x[31],x[32]))
    if str(x[0]) == guide:
        arrayprofileMissmatch = np.array((x[26],x[27],x[28],x[29],x[30],x[31],x[32]),dtype=int)
        
#reading every count file to obtain results for the ecdf and percentile count for annotated type
for line in inExonsCountFile:
    x = line.split('\t')
    exonsMissmatchGlobal.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7]))
    if str(x[0]) == guide:
        arrayexonsMissmatch = np.array((x[1],x[2],x[3],x[4],x[5],x[6],x[7]),dtype=int)

for line in inIntronsCountFile:
    x = line.split('\t')
    intronsMissmatchGlobal.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7]))
    if str(x[0]) == guide:
        arrayintronsMissmatch = np.array((x[1],x[2],x[3],x[4],x[5],x[6],x[7]),dtype=int)

for line in inPromotersCountFile:
    x = line.split('\t')
    promotersMissmatchGlobal.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7]))
    if str(x[0]) == guide:
        arraypromotersMissmatch = np.array((x[1],x[2],x[3],x[4],x[5],x[6],x[7]),dtype=int)

for line in inDNAseCountFile:
    x = line.split('\t')
    DNAseMissmatchGlobal.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7]))
    if str(x[0]) == guide:
        arrayDNAseMissmatch = np.array((x[1],x[2],x[3],x[4],x[5],x[6],x[7]),dtype=int)

for line in inCTCFCountFile:
    x = line.split('\t')
    CTCFMissmatchGlobal.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7]))
    if str(x[0]) == guide:
        arrayCTCFMissmatch = np.array((x[1],x[2],x[3],x[4],x[5],x[6],x[7]),dtype=int)

if(geckoProfile != "no"):
        inGuidesProfile = open(geckoProfile,"r")
        inExonsCountFile = open(geckoExonsCount,"r")
        inIntronsCountFile = open(geckoIntronsCount,"r")
        inPromotersCountFile = open(geckoPromotersCount,"r")
        inDNAseCountFile = open(geckoDNAseCount,"r")
        inCTCFCountFile = open(geckoCTCFCount,"r")
        
        profileMissmatchGlobal.clear()
        exonsMissmatchGlobal.clear()
        intronsMissmatchGlobal.clear()
        promotersMissmatchGlobal.clear()
        DNAseMissmatchGlobal.clear()
        CTCFMissmatchGlobal.clear()

        #reading profile file to obtain results for every mismatch count in the general profile
        next(inGuidesProfile)
        for line in inGuidesProfile:
                x = line.split('\t')
                profileMissmatchGlobal.append((x[26],x[27],x[28],x[29],x[30],x[31],x[32]))

        #reading every count file to obtain results for the ecdf and percentile count for annotated type
        for line in inExonsCountFile:
                x = line.split('\t')
                exonsMissmatchGlobal.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7]))

        for line in inIntronsCountFile:
                x = line.split('\t')
                intronsMissmatchGlobal.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7]))

        for line in inPromotersCountFile:
                x = line.split('\t')
                promotersMissmatchGlobal.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7]))

        for line in inDNAseCountFile:
                x = line.split('\t')
                DNAseMissmatchGlobal.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7]))

        for line in inCTCFCountFile:
                x = line.split('\t')
                CTCFMissmatchGlobal.append((x[1],x[2],x[3],x[4],x[5],x[6],x[7]))

arrayprofileMissmatchGlobal = np.array(profileMissmatchGlobal,dtype=int)
arrayexonsMissmatchGlobal = np.array(exonsMissmatchGlobal,dtype=int)
arrayintronsMissmatchGlobal = np.array(intronsMissmatchGlobal,dtype=int)
arraypromotersMissmatchGlobal = np.array(promotersMissmatchGlobal,dtype=int)
arrayDNAseMissmatchGlobal = np.array(DNAseMissmatchGlobal,dtype=int)
arrayCTCFMissmatchGlobal = np.array(CTCFMissmatchGlobal,dtype=int)


for i in range(0,7):
    #print(str(i)+"MM")
    #ecdf for every annotated type
    ecdfProfile.append(ECDF(arrayprofileMissmatchGlobal[:,i]))
    ecdfExons.append(ECDF(arrayexonsMissmatchGlobal[:,i]))
    ecdfIntrons.append(ECDF(arrayintronsMissmatchGlobal[:,i]))
    ecdfPromoters.append(ECDF(arraypromotersMissmatchGlobal[:,i]))
    ecdfDNAse.append(ECDF(arrayDNAseMissmatchGlobal[:,i]))
    ecdfCTCF.append(ECDF(arrayCTCFMissmatchGlobal[:,i]))
    
    #percentile for every annotated type
    percentileProfile.append(ecdfProfile[i](arrayprofileMissmatch[i]))
    percentileExons.append(ecdfExons[i](arrayexonsMissmatch[i]))
    percentileIntrons.append(ecdfIntrons[i](arrayintronsMissmatch[i]))
    percentilePromoters.append(ecdfPromoters[i](arraypromotersMissmatch[i]))
    percentileDNAse.append(ecdfDNAse[i](arrayDNAseMissmatch[i]))
    percentileCTCF.append(ecdfCTCF[i](arrayCTCFMissmatch[i]))
    #print(ecdfProfile[i](arrayprofileMissmatch[i]))
    #print(ecdfExons[i](arrayexonsMissmatch[i]))
    #print(ecdfIntrons[i](arrayintronsMissmatch[i]))


if len(sys.argv[9])==1:
   # Set data
   df = pd.DataFrame({
   'group': ['A'],
   'EXONS': [percentileExons[missmatch]],
   'GENERAL': [percentileProfile[missmatch]],
   'INTRONS': [percentileIntrons[missmatch]],
   'PROMOTERS': [percentilePromoters[missmatch]],
   'DNAse': [percentileDNAse[missmatch]],
   'CTCF': [percentileCTCF[missmatch]]
   })
   
   # number of variable
   categories=list(df)[1:]
   N = len(categories)

   # We are going to plot the first line of the data frame.
   # But we need to repeat the first value to close the circular graph:
   values=df.loc[0].drop('group').values.flatten().tolist()
   values += values[:1]

   # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
   angles = [n / float(N) * 2 * pi for n in range(N)]
   angles += angles[:1]

   # Initialise the spider plot
   ax = plt.subplot(2,2,1, polar=True)
   #plt.title('RADAR CHART')

   # Draw one axe per variable + add labels labels yet
   plt.xticks(angles[:-1], ['CTCF','DNAse','EXONS','GENERAL','INTRONS','PROMOTERS'], color='black', size=11)

   ax.set_theta_offset(pi / 2)
   ax.set_theta_direction(-1)

   # Draw ylabels
   ax.set_rlabel_position(0)
   plt.yticks([0,0.25,0.50,0.75,1], ["0","0.25","0.50","0.75"], color="grey", size=7)
   plt.ylim(0,1)
   #plt.yticks([10,20,30], ["10","20","30"], color="grey", size=7)
   #plt.ylim(0,40)


   # Plot data
   ax.plot(angles, values, linewidth=1, linestyle='solid')

   # Fill area
   ax.fill(angles, values, 'b', alpha=0.1)

   columns = ('PERCENTILE', 'OFF-TARGETS')
   rows = ('GENERAL', 'EXONS', 'INTRONS','PROMOTERS','DNAse','CTCF')


   offtarget_data=np.vstack((arrayprofileMissmatch[missmatch],arrayexonsMissmatch[missmatch],arrayintronsMissmatch[missmatch],arraypromotersMissmatch[missmatch],arrayDNAseMissmatch[missmatch],arrayCTCFMissmatch[missmatch]))
   percentile_data=np.vstack((float('%.2g' % percentileProfile[missmatch]),float('%.2g' % percentileExons[missmatch]),float('%.2g' % percentileIntrons[missmatch]),float('%.2g' % percentilePromoters[missmatch]),float('%.2g' % percentileDNAse[missmatch]),float('%.2g' % percentileCTCF[missmatch])))
   table_data=np.concatenate((percentile_data,offtarget_data),axis=1)

   plt.subplot(2, 2, 2)
   #plt.title('SUMMARY TABLE')
   table=plt.table(cellText=table_data,rowLabels=rows,colLabels=columns,loc='center',colWidths=[0.35 for x in columns])
   #table=plt.table(cellText=table_data,rowLabels=rows,colLabels=columns,loc='center',bbox=[0,-0.65,1,0.65])
   #table.set_fontsize(20)
   table.auto_set_font_size(False)
   table.set_fontsize(10)
   table.scale(1, 1.5)
   #table.scale(1.5, 1.5)
   #cellDict=table.get_celld()
   #cellDict[(0,0)].set_width(0.1)plt-
   plt.axis('off')


   #datacount=(446,525,356,710,478,2336,1319,532,449,589,2016,383,771,2573,576,439,644,366,2226,2451)
   datacount=arrayguidesExtendedProfile[missmatch*5]/(max(arrayguidesExtendedProfile[missmatch*5]))
   data = np.array(datacount,dtype=float)
   data = np.around(data,decimals=3)
   data.shape=(1,len(datacount))

   string = guide[0:20]


   A = arrayguidesExtendedProfile[missmatch*5+1]/(max(arrayguidesExtendedProfile[missmatch*5]))
   C = arrayguidesExtendedProfile[missmatch*5+2]/(max(arrayguidesExtendedProfile[missmatch*5]))
   G = arrayguidesExtendedProfile[missmatch*5+3]/(max(arrayguidesExtendedProfile[missmatch*5]))
   T = arrayguidesExtendedProfile[missmatch*5+4]/(max(arrayguidesExtendedProfile[missmatch*5]))

   ind = np.arange(0,len(string),1) + 0.15  # the x locations for the groups
   width = 0.7 # the width of the bars: can also be len(x) sequence

   motif=plt.subplot(2,1,2)
   p1 = plt.bar(ind, A,width,color='#d62728',align='edge')
   p2 = plt.bar(ind, C,width,bottom=A,align='edge')
   p3 = plt.bar(ind, G,width,bottom=A+C,align='edge')
   p4 = plt.bar(ind, T,width,bottom=C+G+A,align='edge')

   plt.xlim(0,len(string))
   plt.xticks([])
   #plt.ylabel('Missmatches')
   #plt.title('MOTIF LOGO')

   plt.legend((p1[0], p2[0],p3[0],p4[0]), ('A', 'C','G', 'T'))


   #plt.subplot(2, 1, 2)
   #table=plt.table(cellText=data,rowLabels=row_string,colLabels=string,loc='bottom',colWidths=[0.05 for x in string])
   table=plt.table(cellText=data,colLabels=string,loc='bottom')
   #plt.xticks(ind)
   #table.set_fontsize(20)
   table.auto_set_font_size(False)
   table.set_fontsize(10)

   #table.scale(1.5, 1.5)
   #cellDict=table.get_celld()
   #cellDict[(0,0)].set_width(0.1)
   #plt.axis('off')

   #plt.tight_layout()

   #plt.suptitle(string +' '+ "ANALYSIS",horizontalalignment='center', color='black', weight='bold',size='large')
   plt.suptitle(str(missmatch)+"MM",horizontalalignment='center', color='black', weight='bold',size='large')

   plt.show()
   # generate your plot
   #plt.savefig("myfig.png",dpi=400,bbox_inches='tight')

   #plt.savefig("image.png",bbox_inches='tight',dpi=400)

else:

   def make_spider(row, title,count):

      df = pd.DataFrame({
      'group': ['A'],
      'EXONS': [percentileExons[row]],
      'GENERAL': [percentileProfile[row]],
      'INTRONS': [percentileIntrons[row]],
      'PROMOTERS': [percentilePromoters[row]],
      'DNAse': [percentileDNAse[row]],
      'CTCF': [percentileCTCF[row]]
      })

 
    # number of variable
      categories=list(df)[1:]
      N = len(categories)
      
      # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
      angles = [n / float(N) * 2 * pi for n in range(N)]
      angles += angles[:1]
      
      # Initialise the spider plot
      ax = plt.subplot(2,(uppermm-lowermm)+1,count, polar=True)
      #table=plt.table(cellText=data,colLabels=string,loc='bottom')
      
      # If you want the first axis to be on top:
      ax.set_theta_offset(pi / 2)
      ax.set_theta_direction(-1)
      
      # Draw one axe per variable + add labels labels yet
      plt.xticks(angles[:-1], ['CTCF'+' ('+str(arrayCTCFMissmatch[row])+')'+' ('+str(float('%.2g' % percentileCTCF[row]))+')','DNAse'+' ('+str(arrayDNAseMissmatch[row])+')'+' ('+str(float('%.2g' % percentileDNAse[row]))+')','EXONS'+' ('+str(arrayexonsMissmatch[row])+')'+' ('+str(float('%.2g' % percentileExons[row]))+')','GENERAL'+' ('+str(arrayprofileMissmatch[row])+')'+' ('+str(float('%.2g' % percentileProfile[row]))+')','INTRONS'+' ('+str(arrayintronsMissmatch[row])+')'+' ('+str(float('%.2g' % percentileIntrons[row]))+')','PROMOTERS'+' ('+str(arraypromotersMissmatch[row])+')'+' ('+str(float('%.2g' % percentilePromoters[row]))+')'], color='black', size=10)
      
      # Draw ylabels
      ax.set_rlabel_position(0)
      plt.yticks([0,0.25,0.50,0.75], ["0","0.25","0.50","0.75"], color="black", size=9)
      plt.ylim(0,1)
      
      # Ind1
      values=df.loc[0].drop('group').values.flatten().tolist()
      values += values[:1]
      ax.plot(angles, values, linewidth=2, linestyle='solid')
      ax.fill(angles, values, alpha=0.4)

      colors = ["white", "white","white", "white","white","white"]
      texts = ["CTCF:"+' '+str(max(arrayCTCFMissmatchGlobal[:,row])),"DNAse:"+' '+str(max(arrayDNAseMissmatchGlobal[:,row])),"EXONS:"+' '+str(max(arrayexonsMissmatchGlobal[:,row])),"GENERAL:"+' '+str(max(arrayprofileMissmatchGlobal[:,row])),"INTRONS:"+' '+str(max(arrayintronsMissmatchGlobal[:,row])),"PROMOTERS:"+' '+str(max(arraypromotersMissmatchGlobal[:,row]))]
      patches = [mpatches.Patch(color=colors[i] ,label="{:s}".format(texts[i]) ) for i in range(len(texts))]
      plt.legend(handles=patches, loc=(-0.35,0.85),labelspacing=0.5, fontsize='10', ncol=1,handlelength=0, handletextpad=0,title="MAX VALUE")     

#       plt.gcf().text(0.05+((count-1)/(1+uppermm-lowermm)), 0.90, "CTCF"+str(max(arrayCTCFMissmatchGlobal[:,row])), fontsize=10)
#       plt.gcf().text(0.05+((count-1)/(1+uppermm-lowermm)), 0.87, "CTCF"+str(max(arrayCTCFMissmatchGlobal[:,row])), fontsize=10)
#       plt.gcf().text(0.05+((count-1)/(1+uppermm-lowermm)), 0.84, "CTCF"+str(max(arrayCTCFMissmatchGlobal[:,row])), fontsize=10)
#       plt.gcf().text(0.05+((count-1)/(1+uppermm-lowermm)), 0.81, "CTCF"+str(max(arrayCTCFMissmatchGlobal[:,row])), fontsize=10)
#       plt.gcf().text(0.05+((count-1)/(1+uppermm-lowermm)), 0.78, "CTCF"+str(max(arrayCTCFMissmatchGlobal[:,row])), fontsize=10)

      # Add a title
      plt.title(title, size=8, y=1.1)


   def make_motif(row,count):

      datacount=arrayguidesExtendedProfile[row*5]/(max(arrayguidesExtendedProfile[row*5]))
      data = np.array(datacount,dtype=float)
      data = np.around(data,decimals=3)
      data.shape=(1,len(datacount))

      string = guide[0:20]

      A = arrayguidesExtendedProfile[row*5+1]/(max(arrayguidesExtendedProfile[row*5]))
      C = arrayguidesExtendedProfile[row*5+2]/(max(arrayguidesExtendedProfile[row*5]))
      G = arrayguidesExtendedProfile[row*5+3]/(max(arrayguidesExtendedProfile[row*5]))
      T = arrayguidesExtendedProfile[row*5+4]/(max(arrayguidesExtendedProfile[row*5]))


      ind = np.arange(0,len(string),1) + 0.15  # the x locations for the groups
      width = 0.7 # the width of the bars: can also be len(x) sequence

      # motif=plt.subplot(3,4,row)
      # p1 = plt.bar(ind, A,width,color='#d62728',align='edge')
      # p2 = plt.bar(ind, C,width,bottom=A,align='edge')
      # p3 = plt.bar(ind, G,width,bottom=A+C,align='edge')
      # p4 = plt.bar(ind, T,width,bottom=C+G+A,align='edge')

      
      # #plt.ylabel('Missmatches')
      # #plt.title('MOTIF LOGO')

      # plt.legend((p1[0], p2[0],p3[0],p4[0]), ('A', 'C','G', 'T'))

      plt.subplot(2,(uppermm-lowermm)+1,count+(uppermm-lowermm)+1)

      p1 = plt.bar(ind, A,width,color='#d62728',align='edge')
      p2 = plt.bar(ind, C,width,bottom=A,align='edge')
      p3 = plt.bar(ind, G,width,bottom=A+C,align='edge')
      p4 = plt.bar(ind, T,width,bottom=C+G+A,align='edge')

      plt.legend((p1[0], p2[0],p3[0],p4[0]), ('A', 'C','G', 'T'))

      plt.xlim(0,len(string))
      plt.xticks([])
      table=plt.table(cellText=data,colLabels=string,loc='bottom')
      #plt.xticks(ind)
      #table.set_fontsize(20)
      table.auto_set_font_size(False)
      table.set_fontsize(4.5)

   # ------- PART 2: Apply to all individuals
   # initialize the figure
   my_dpi=96
   plt.figure(figsize=(1000/my_dpi, 1000/my_dpi), dpi=my_dpi)

   # Create a color palette:
   #my_palette = plt.cm.get_cmap("Set2",(uppermm-lowermm)+1)

   count=1
   # Loop to plot
   for row in range(lowermm, uppermm+1):
      make_spider(row=row, title=str(row) +' MISMATCHES',count=count)
      make_motif(row=row,count=count)
      #plt.subplot(3,4,count*2)
      count=count+1
      
   #plt.tight_layout()
   plt.subplots_adjust(top=0.9,bottom=0.05,left=0.05,right=0.95,wspace=0.1)
   plt.show()

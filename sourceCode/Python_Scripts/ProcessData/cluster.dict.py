#!/usr/bin/env python

# Uso tuple invece di liste, joino le parti delle lines che non uso, cambio str in int 
# Input 2.8 -> 17 Gb, doppio del tempo per clusterizzare, ma sempre sotto i 10 minuti
# Input di 5.6 gb -> non funziona
#sys1 è target file
#sys2 is 'addGuide' or 'no' -> only for web server, only for search with only ref
#sys3 is True to keep column 5 (Pos cluster) and 9 (Total) and added guide, False to do clusterization but do not report the added columns
#sys4 is True if cluster only (no append of Total column or adding of Cluster position, because already present), False otherwise
#sys5 is guides.txt for slow clustering
#sys6 is result name
#sys [-1] can be result name OR 'Step [1/5]', needed to print current operation
#So if 'Step' in sys.argv[-1] do a print of current operation
#if sys contains 'orderChr', clusters will be ordered by chr
#Output column (not written): Bulge_type, Guide, Target, chr, pos, pos_cluster (optional), direction, mms, bulge, total(optional), real guide(optional)

#If sys contains orderChr, first call subprocess to sort by chr, cluster pos, real guide
#Poi metto in una lista i target che hanno lo stesso chr clusterpos direction real guide, faccio sort per total e mms, poi lo salvo e passo al prossimo
#Ottengo così i cluster ordinati per chr, ma non ordinati tra di loro (eg un cluster con total 5 può apparire prima di uno con tot

#NOTE with new search, cluster position is already present, so code for that column is commented
#NOTE 06/03 PAM  -> removed PAM Disruption calculation
#NOTE 11/03 Removed TOTAL column (already present in search phase)
import time
import sys
import subprocess
import os 
MAX_LIMIT = 5000000000
alphabet = '-RYSWKMBDHVryswkmbdhvACGTacgt'      #For forcing IUPAC before reference target
pam_at_end = True
guides_to_check = set()        #Set of error guides, do not do computation on them
if os.path.exists('./guides_error.txt'):
    with open('guides_error.txt') as g_e:
        for line in g_e:
            line = line.strip()
            guides_to_check.add(line)
            if line[0] == 'N':
                pam_at_end = False
with open(sys.argv[5]) as guides:
    for guide in guides:
        if guide[0] == 'N':
            pam_at_end = False

strand_to_check = '-'
if pam_at_end:
    strand_to_check = '+'

start = time.time()
total_targets = []
guides_dict = dict()
addGuide = False                #Add real guide to last column for better grep
if 'addGuide' in sys.argv[:]:
    addGuide = True
if sys.argv[3] == 'True':
    keep_columns = True
else:
    keep_columns = False

add_for_final = '\tSamples\tAnnotation Type\tReal Guide\n'

result_name = sys.argv[6] + '.cluster.txt'
cluster_only = False
if sys.argv[4] == 'True':
    cluster_only = True

process = subprocess.Popen(['wc', '-l', sys.argv[1]], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = process.communicate()
total_line = int(out.decode('UTF-8').split(' ')[0])
step = sys.argv[-1]
#Do cluster by orderChr NEW VERSION
if 'orderChr' in sys.argv[:]:
    if 'Step' in step:
        print(step + ': Sorting Targets', end = '\r')
    subprocess.run(['sort -k4,4 -k6,6 -k15,15 -k7,7 ' + sys.argv[1] + ' > ' + result_name + '.tmp_sort.txt'], shell = True)      #sort input file by chr, clusterpos, real guide, direction
    
    start_time = time.time()
    if 'Step' in step:
        print(step + ': Clustering Targets', end = '\r')
    with open (result_name + '.tmp_sort.txt') as targets, open(result_name, 'w+') as result:
        #Write Header    
        if 'total' not in sys.argv[:]:    
            if addGuide:
                if keep_columns:
                    result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tReal Guide\n')
                else:
                    result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge_Size\tReal Guide\n')
            else:
                if keep_columns:
                    result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\n')
                else:
                    result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge_Size\n')
        else:
            result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq' + add_for_final)
        current_chr_pos_dir_g = ''
        current_cluster = []
        for line in targets:
            if '#' in line:     #Skip header
                continue
            line = line.strip().split('\t')
            if line[3] + line[5] + line[6] + line[14] != current_chr_pos_dir_g:      #NEW CLUSTER FOUND
                #Sort previous cluster
                current_cluster.sort(key = lambda x: (x[9], x[7], [alphabet.index(c) for c in x[2]]))    #Sort by total, mismatches, and prioritize IUPAC characters
                
                #Save previous cluster
                for t in current_cluster:
                    result.write('\t'.join(t) + '\n')

                #Initialize new cluster
                current_cluster = []
                current_cluster.append(line)
                current_chr_pos_dir_g = line[3] + line[5] + line[6] + line[14]
            else:                                                       #IN THE SAME CLUSTER
                current_cluster.append(line)
        #Sort and save last cluster
        current_cluster.sort(key = lambda x: (x[9], x[7], [alphabet.index(c) for c in x[2]]))
        for t in current_cluster:
            result.write('\t'.join(t) + '\n')
        
        #Remove tmp sort
        os.remove(result_name + '.tmp_sort.txt')
        # print('END Clustering (Local)', time.time() - start_time)
        exit()

if total_line > MAX_LIMIT:
    print('Max limit file reached, using clustering by single guide')
    ok_guides = []
    write_header = True
    with open(sys.argv[5]) as guides:
        for guide in guides:
            if guide in guides_to_check:
                continue
            current_count = 0
            cluster_ok = True
            guide = guide.strip()
            #DO SLOW CLUSTERING       
            with open (sys.argv[1]) as targets:
                for line in targets:
                    line = line.strip().split('\t')
                    if '#' in line[0] or line[1].replace('-','') != guide:
                        continue
                    # if not cluster_only:
                    #     line.append(str(int(line[7]) + int(line[8])))   #Total column
                        
                    current_count += 1
                    # if current_count > MAX_LIMIT:
                    #     print('The guide ' + guide + ' has more than ' + str(MAX_LIMIT) + ' targets. Skipping...')
                    #     with open('./guides_error.txt', 'a+') as guides_error:
                    #         guides_error.write(guide + '\n')
                    #     cluster_ok = False
                    #     del guides_dict[guide]
                    #     break
                    try:
                        guides_dict[line[1].replace('-','')].append(('\t'.join(line[:2]), line[2], line[3], int(line[4]), int(line[5]), str(line[6]), int(line[7]), int(line[8]), int(line[9]), '\t'.join(line[10:])))    #[('type\tguide', 'target', 'chr', pos, clusterpos, 'dir', mm, bul, tot)]
                    except:
                        guides_dict[line[1].replace('-','')] = [('\t'.join(line[:2]), line[2], line[3], int(line[4]), int(line[5]), str(line[6]), int(line[7]), int(line[8]), int(line[9]),'\t'.join(line[10:]) )]
            if not cluster_ok:
                continue
            # if not cluster_only:       
            #     # print('Created \'Total\' and \'Position Cluster\' columns:', time.time() - start)
            # else:
            #     # print('Loaded targets: ', time.time() - start)
            start = time.time()
            # total_targets.sort(key = lambda x: ( x[3] , int(x[-1]) ))
            for k in guides_dict.keys():
                guides_dict[k].sort(key = lambda x: ( x[2] , x[4] ))    #Order by chr_clusterpos
            #total_targets.sort(key = lambda x: ( x[3] , int(x[5]) ))

            # print('Targets sorted:', time.time() - start)

            # print('Start clustering')
            start_time = time.time()

            with open(result_name, 'a+') as result:
                if write_header:    
                    if 'total' not in sys.argv[:]:      #TODO fix for offline release, praticamente se sto facendo cluster su total.txt metto l'header custom
                        if addGuide:
                            if keep_columns:
                                result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tReal Guide\n')
                            else:
                                result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge_Size\tReal Guide\n')
                        else:
                            if keep_columns:
                                result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\n')
                            else:
                                result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge_Size\n')
                    else:
                        result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq' + add_for_final)
                    write_header = False
                total_targets = []
                for k in guides_dict.keys():
                    total_targets += guides_dict[k]
                total_list = []

                first_line = total_targets[0]
                # current_chr_pos = first_line[3] + ' ' + first_line[9]
                current_chr_pos = first_line[2] + str(first_line[4]) + first_line[5]  #chr clusterpos direction

                total_list.append([first_line])

                for line in total_targets[1:]:
                    #if line[3] + ' ' + line[9] != current_chr_pos:
                    if line[2] + str(line[4]) + line[5] != current_chr_pos:
                        # total_list[-1].sort(key = lambda x: int(x[8]))
                        total_list[-1].sort(key = lambda x: (x[8], x[6], [alphabet.index(c) for c in x[1]]))   #Order cluster by total and mms, and prioritize targets with IUPAC 
                        total_list.append([line])
                        # current_chr_pos = line[3] + ' ' + line[9]
                        current_chr_pos = line[2] + str(line[4]) + line[5]
                    else:
                        total_list[-1].append(line)     

                total_list[-1].sort(key = lambda x: (x[8], x[6], [alphabet.index(c) for c in x[1]]))       #Order last cluster by total and mms, and prioritize targets with IUPAC

                total_list.sort(key = lambda x: x[0][8])        #Order all clusters by total of top1
                if 'orderChr' in sys.argv[:]:
                    total_list.sort(key = lambda x: x[0][2])    # and then for chr, needed for sample annotation and summary by position
                if addGuide:
                    if keep_columns:
                        for cluster in total_list:
                            for target in cluster:
                                result.write('\t'.join(str(x) for x in target).strip() + '\t' + target[0].split('\t')[1].replace('-','') + '\n')
                    else:
                        for cluster in total_list:
                            for target in cluster:
                                result.write('\t'.join(str(x) for x in target[:4]) + '\t' + '\t'.join(str(x) for x in target[5:8]) +'\t' + target[9].strip() + '\t' + target[0].split('\t')[1].replace('-','') + '\n')

                else:
                    if keep_columns:
                        for cluster in total_list:
                            for target in cluster:
                                result.write('\t'.join([str(x) for x in target]).strip() + '\n')
                    else:
                        for cluster in total_list:
                            for target in cluster:
                                result.write('\t'.join(str(x) for x in target[:4]) + '\t' + '\t'.join(str(x) for x in target[5:8]) +'\t' + target[9].strip() + '\n')
            del guides_dict[guide]
            ok_guides.append(guide)
            # print("Clustering runtime: %s seconds" % (time.time() - start_time))
    with open(sys.argv[5], 'w') as guides:
        guides.write('\n'.join(ok_guides))        

else:
    #####DO FAST CLUSTERING
    with open (sys.argv[1]) as targets:
        for line in targets:
            if '#' in line:
                continue
            line = line.strip().split('\t')
            if line[1].replace('-','') in guides_to_check:
                continue
            # if not cluster_only:
            #     line.append(str(int(line[7]) + int(line[8])))   #Add Total column
                
            try:
                guides_dict[line[1].replace('-','')].append(('\t'.join(line[:2]), line[2], line[3], int(line[4]), int(line[5]), str(line[6]), int(line[7]), int(line[8]), int(line[9]), '\t'.join(line[10:])))    #[('type\tguide', 'target', 'chr', pos, clusterpos, 'dir', mm, bul, tot)]
            except:
                guides_dict[line[1].replace('-','')] = [('\t'.join(line[:2]), line[2], line[3], int(line[4]), int(line[5]), str(line[6]), int(line[7]), int(line[8]), int(line[9]),'\t'.join(line[10:]) )]
            #total_targets.append(line)
    # if not cluster_only:       
    #     # print('Created \'Total\' and \'Position Cluster\' columns:', time.time() - start)
    # else:
    #     # print('Loaded targets: ', time.time() - start)
    start = time.time()
    # total_targets.sort(key = lambda x: ( x[3] , int(x[-1]) ))
    for k in guides_dict.keys():
        guides_dict[k].sort(key = lambda x: ( x[2] , x[4] ))        #Order by chr_clusterpos
    #total_targets.sort(key = lambda x: ( x[3] , int(x[5]) ))

    # print('Targets sorted:', time.time() - start)

    # print('Start clustering')
    start_time = time.time()

    with open(result_name, 'w+') as result:
        if 'total' not in sys.argv[:]:      #TODO fix for offline release, praticamente se sto facendo cluster su total.txt metto l'header custom
            if addGuide:
                if keep_columns:
                    result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tReal Guide\n')
                else:
                    result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge_Size\tReal Guide\n')
            else:
                if keep_columns:
                    result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\n')
                else:
                    result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge_Size\n')
        else:
            result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq' + add_for_final)

        total_targets = []
        for k in guides_dict.keys():
            total_targets += guides_dict[k]
        total_list = []

        if not total_targets:
            print('No targets to clusterize, exit...')
            sys.exit()
        first_line = total_targets[0]
        # current_chr_pos = first_line[3] + ' ' + first_line[9]
        current_chr_pos = first_line[2] + str(first_line[4]) + first_line[5]

        total_list.append([first_line])

        for line in total_targets[1:]:
            #if line[3] + ' ' + line[9] != current_chr_pos:
            if line[2] + str(line[4]) + line[5] != current_chr_pos:
                # total_list[-1].sort(key = lambda x: int(x[8]))
                total_list[-1].sort(key = lambda x: (x[8], x[6], [alphabet.index(c) for c in x[1]]))   #Order cluster by total and mms 
                total_list.append([line])
                # current_chr_pos = line[3] + ' ' + line[9]
                current_chr_pos = line[2] + str(line[4]) + line[5]
            else:
                total_list[-1].append(line)     

        total_list[-1].sort(key = lambda x: (x[8], x[6], [alphabet.index(c) for c in x[1]]))       #Order last cluster by total and mms 
        total_list.sort(key = lambda x: x[0][8])        #Order all clusters by total of top1
        if 'orderChr' in sys.argv[:]:
            total_list.sort(key = lambda x: x[0][2])    # and then for chr, needed for sample annotation and summary by position
        if addGuide:
            if keep_columns:
                for cluster in total_list:
                    for target in cluster:
                        result.write('\t'.join(str(x) for x in target).strip() + '\t' + target[0].split('\t')[1].replace('-','') + '\n')
            else:
                for cluster in total_list:
                    for target in cluster:
                        result.write('\t'.join(str(x) for x in target[:4]) + '\t' + '\t'.join(str(x) for x in target[5:8]) +'\t' + target[9].strip() + '\t' + target[0].split('\t')[1].replace('-','') + '\n')

        else:
            if keep_columns:
                for cluster in total_list:
                    for target in cluster:
                        result.write('\t'.join([str(x) for x in target]).strip() + '\n')
            else:
                for cluster in total_list:
                    for target in cluster:
                        result.write('\t'.join(str(x) for x in target[:4]) + '\t' + '\t'.join(str(x) for x in target[5:8]) +'\t' + target[9].strip() + '\n')

    # print("Clustering runtime: %s seconds" % (time.time() - start_time))
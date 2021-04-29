#!/usr/bin/env python

# Python Program to search CRISPR/Cas complex into a genome
import subprocess					# run c++ executable
import time
import os							# instructions manage directories
from os import listdir
from os.path import isfile, join, isdir
import shutil						# remove directory and its content
import sys							# input argv
from subprocess import Popen, PIPE
import glob
import pandas as pd
import multiprocessing
# path where this file is located
origin_path = os.path.dirname(os.path.realpath(__file__))
# conda path
conda_path = "opt/crispritz/"

if '--debug' in sys.argv[1:]:
    #for quick local tests
    corrected_origin_path = origin_path+'/sourceCode/'
else:
     # path corrected to use with conda
    corrected_origin_path = origin_path[:-3]+conda_path


def checkExistance(f_path, element):
    if element == 'f':  # check file
        if not isfile(f_path):
            print('ERROR! ' + f_path + ' is not a file.')
            sys.exit()
    else:  # check directory
        if not isdir(f_path):
            print('ERROR! ' + f_path + ' is not a directory.')
            sys.exit()

# ALL ACTIVE FUNCTIONS IN CRISPRITZ


def indexGenome():
    if (len(sys.argv) < 7 or 'help' in sys.argv[1:]):
        print("WARNING: Too few arguments to function index-genome. Please provide:\n",
              "\nEXAMPLE CALL: crispritz.py index-genome name_genome genomeDirectory(FASTA)/ pamFile.txt -bMax 2\n",
              "\n<name_genome>: Name of the genome to create",
              "\n<genomeDirectory>: Directory containing a genome in .fa or .fasta format, need to be separated into single chromosome files.",
              "\n<pamFile>: Text file containing the PAM (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence",
              "\n-bMax <maxBulges_num>: Number of bulges allowed for the search phase",
              "\n-th <num_thread>: (Optional) Number of threads to use. Default uses all of the available threads")
        sys.exit()

    nameGenome = sys.argv[2]								# save name of the genome
    dirGenome = os.path.realpath(sys.argv[3])
    checkExistance(dirGenome, 'd')
    dirPAM = os.path.realpath(sys.argv[4])
    checkExistance(dirPAM, 'f')
    filePAM = open(os.path.realpath(sys.argv[4]), "r")
    listChrs = os.listdir(dirGenome)						# save list of chromosomes

    max_bulges = (sys.argv).index("-bMax") + 1
    try:
        max_bulges = sys.argv[max_bulges]
    except IndexError:
        print('ERROR! No value provided for the -bMax option, please provide a number')
        sys.exit()
    try:
        int(max_bulges)
    except ValueError:
        print('ERROR! The value for the -bMax option is invalid, please provide a number')
        sys.exit()
    # retrive PAM
    PAM = filePAM.read()
    PAM_size = int(PAM.split()[1])
    if(PAM_size < 0):
        PAM_size = PAM_size * -1
        PAM = PAM.split()[0][0:PAM_size]
    else:
        PAM = PAM.split()[0][-PAM_size:]

    TSTgenome = PAM + "_" + max_bulges + "_" + \
        nameGenome					# name of the genome in TST format
    dirTSTgenome = "./genome_library/" + TSTgenome		# dir of the genome in TST format

    print(TSTgenome, "Indexing generation:")

    # read number of threads
    th = 1000
    if "-th" in sys.argv[1:]:
        try:
            th = (sys.argv).index("-th") + 1
            th = int(sys.argv[th])
        except:
            print(
                "ATTENTION! Check the thread option: -th <th_num> (th_num is an integer)")
            sys.exit()

    if os.path.isdir(dirTSTgenome):							# check if TSTgenome dir exists
        shutil.rmtree(dirTSTgenome)							# remove old TSTgenome dir
    os.makedirs(dirTSTgenome)								# build new TSTgenome dir
    os.chdir(dirTSTgenome)									# move into the TSTgenome dir

    # run buildTST
    start_time = time.time()
    for f in listChrs:
        if f.strip().split('.')[-1] != 'fa':
            continue
        print("Indexing:", f)
        subprocess.run([corrected_origin_path+"buildTST",
                        str(dirGenome)+"/"+str(f), str(dirPAM), str(th), max_bulges])
    print("Finish indexing")
    print("Indexing runtime: %s seconds" % (time.time() - start_time))


def searchTST():
    if (len(sys.argv) < 10 or 'help' in sys.argv[1:]):
        print("WARNING: Too few arguments to function search with bulges selected\n",
              "\ncrispritz.py search help to see paramaters and example call\n")
        sys.exit()

    nameGenome = os.path.realpath(sys.argv[2])
    PAM = os.path.realpath(sys.argv[3])								# save PAM
    fileGuide = os.path.realpath(sys.argv[4])
    nameResult = (sys.argv[5])  # name of result file
    dirTSTgenome = os.path.realpath(sys.argv[2])+"/"
    try:
        max_bulges = dirTSTgenome.split("/")[-2].split("_")[1]
    except IndexError:
        print('ERROR! The directory name must contain the PAM and the maximum available bulges, please create the Genome Index using crispritz')
        sys.exit()

    checkExistance(PAM, 'f')
    checkExistance(fileGuide, 'f')
    if not os.path.isdir(dirTSTgenome):				# check if TSTgenome dir exists
        print("ATTENTION! You have to generate the index of \"" + nameGenome +
              "\" with \"" + PAM + "\", before the search using index!")
        sys.exit()

    # read number of mismatches
    mm = 0
    if "-mm" in sys.argv[1:]:
        try:
            mm = (sys.argv).index("-mm") + 1
            mm = int(sys.argv[mm])
        except:
            print(
                "ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is a number)")
            sys.exit()

    # read number of bulge RNA
    bRNA = 0
    if "-bRNA" in sys.argv[1:]:
        try:
            bRNA = (sys.argv).index("-bRNA") + 1
            bRNA = int(sys.argv[bRNA])
        except:
            print(
                "ATTENTION! Check the bulge RNA option: -bRNA <bRNA_num> (bRNA_num is a number)")
            sys.exit()

    # read number of bulge DNA
    bDNA = 0
    if "-bDNA" in sys.argv[1:]:
        try:
            bDNA = (sys.argv).index("-bDNA") + 1
            bDNA = int(sys.argv[bDNA])
        except:
            print(
                "ATTENTION! Check the bulge DNA option: -bDNA <bDNA_num> (bDNA_num is a number)")
            sys.exit()

    if (int(max_bulges) < bDNA):
        print("WARNING! Max available bulges (" + max_bulges +
              ") is smaller than input DNA bulges (" + str(bDNA) + "). Using " + max_bulges + " bulges")
        bDNA = max_bulges
    if (int(max_bulges) < bRNA):
        print("WARNING! Max available bulges (" + max_bulges +
              ") is smaller than input RNA bulges (" + str(bRNA) + "). Using " + max_bulges + " bulges")
        bRNA = max_bulges
    # read number of threads
    th = 1
    if "-th" in sys.argv[1:]:
        try:
            th = (sys.argv).index("-th") + 1
            th = int(sys.argv[th])
        except:
            print(
                "ATTENTION! Check the thread option: -th <th_num> (th_num is an integer)")
            sys.exit()

    # results writing
    r = "no"
    if "-r" in sys.argv[1:]:
        r = "r"
    if "-p" in sys.argv[1:]:
        r = "p"
    if "-t" in sys.argv[1:]:
        r = "t"
    if r == "no":
        print("Please select an output")
        sys.exit()

    # Check '-scores' directory input
    if "-scores" in sys.argv[1:]:
        try:
            idx_genome_fasta = (sys.argv).index("-scores") + 1
            idx_genome_fasta = os.path.realpath(sys.argv[idx_genome_fasta])
        except:
            print(
                "ERROR! Please select the directory containing the fasta files of the genome")
            sys.exit()
        checkExistance(idx_genome_fasta, 'd')
    # Check input correctness
    file_correct_ext = [f for f in listdir(dirTSTgenome) if isfile(join(
        dirTSTgenome, f)) and not f.endswith('.bin')]  # Get files not ending with .bin
    if len(file_correct_ext) != 0:		# Some files do not have .bin
        print('ERROR! The directory contains files that are not compatible with the selected search type. Please ensure that the input directory contains only .bin files')
        sys.exit()

    # run searchOnTST
    print("Search START")
    start_time = time.time()
    subprocess.run([corrected_origin_path+"searchTST", str(dirTSTgenome), str(fileGuide),
                    str(mm), str(bDNA), str(bRNA), str(PAM), str(nameResult), str(r), str(th), max_bulges])
    print("Search END")
    print("Search runtime: %s seconds" % (time.time() - start_time))
    if "-scores" in sys.argv[1:]:

        try:
            idx_genome_fasta = (sys.argv).index("-scores") + 1
            idx_genome_fasta = os.path.realpath(sys.argv[idx_genome_fasta])
        except:
            print(
                "ERROR! Please select the directory containing the fasta files of the genome")
            sys.exit()

        pam_guide = len(open(PAM).readline().split(" ")[0])
        pam_at_beginning = int(open(PAM).readline().split(" ")[1])
        if (pam_guide != 23 or pam_at_beginning < 0):  # Also block scoring pam at beginning
            print("WARNING: The CFD score and the Doench score can be calculated only for guides with 20bp and a 3bp PAM (SpCas9)")

        target_filename = os.path.realpath(nameResult)
        subprocess.run([corrected_origin_path+'Python_Scripts/Scores/scores.py',
                        target_filename + '.targets.txt', idx_genome_fasta + "/", str(PAM), str(fileGuide)])


def searchBruteForce():
    if (len(sys.argv) < 9 or 'help' in sys.argv[1:]):
        print("WARNING: Too few arguments to function search. Please provide:\n",
              "\nEXAMPLE CALL MISMATCHES ONLY: crispritz.py search genomeDirectory/ pamFile.txt guidesFile.txt outputFile -mm 4 -th 4 -scores genomeDirectory(FASTA)/ -t\n",
              "\nEXAMPLE CALL MISMATCHES + BULGES: crispritz.py search indexGenomeDirectory/ pamFile.txt guidesFile.txt outputFile -mm 4 -bDNA 1 -bRNA 2 -th 4 -scores genomeDirectory(FASTA)/ -t\n",
              "\n<genomeDirectory>: Directory containing a genome in .fa or .fasta format (.bin format if bulges present), need to be separated into single chromosome files",
              "\n<pamFile>: Text file containing the PAM sequence (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence",
              "\n<guidesFile>: Text file containing one or more guides (including a number of Ns equal to the length of the PAM sequence)",
              "\n<outputFile>: Name of output file",
              "\n-mm <mm_num>: Number of allowed mismatches",
              "\n-bRNA <bRNA_num>: (Optional) Size of RNA bulges",
              "\n-bDNA <bDNA_num>: (Optional) Size of DNA bulges",
              "\n-th < num_thread >: (Optional) Number of threads to use. Default uses all of the available threads (ONE for bulge search)",
              "\n-scores <genomeDirectoryInFastaFormat>: (Optional) Directory containing the genome in .fa or .fasta format, necessary to extract sequences for Doench Score Function",
              "\n{-r,-p,-t}: Output type (-r off-targets list only, -p profile only, -t off-targets AND profile)\n")
        sys.exit()

    genomeDir = os.path.realpath(sys.argv[2])+"/"
    filePAM = os.path.realpath(sys.argv[3])
    fileGuide = os.path.realpath(sys.argv[4])
    result = sys.argv[5]

    # read number of mismatches
    mm = 0
    if "-mm" in sys.argv[1:]:
        try:
            mm = (sys.argv).index("-mm") + 1
            mm = int(sys.argv[mm])
        except:
            print(
                "ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is an integer)")
            sys.exit()

    # read number of mismatches
    th = 1000
    if "-th" in sys.argv[1:]:
        try:
            th = (sys.argv).index("-th") + 1
            th = int(sys.argv[th])
        except:
            print(
                "ATTENTION! Check the thread option: -th <th_num> (th_num is an integer)")
            sys.exit()

    # results writing
    r = "no"
    if "-r" in sys.argv[1:]:
        r = "r"
    if "-p" in sys.argv[1:]:
        r = "p"
    if "-t" in sys.argv[1:]:
        r = "t"
    if r == "no":
        print("Please select an output")
        sys.exit()

    # variant
    variant = 1  # always search with IUPAC variant

    # Check '-scores' directory
    if "-scores" in sys.argv[1:]:
        try:
            idx_genome_fasta = (sys.argv).index("-scores") + 1
            idx_genome_fasta = os.path.realpath(sys.argv[idx_genome_fasta])
        except:
            print(
                "ERROR! Please select the directory containing the fasta files of the genome")
            sys.exit()
        checkExistance(idx_genome_fasta, 'd')

    # Check input correctness
    checkExistance(genomeDir, 'd')
    checkExistance(filePAM, 'f')
    checkExistance(fileGuide, 'f')

    file_correct_ext = [f for f in listdir(genomeDir) if isfile(join(genomeDir, f)) and not (
        f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.fai'))]

    if len(file_correct_ext) != 0:  # Some file other than .fa, .fasta or .fai are present in the directory
        print('ERROR! The directory contains files that are not compatible with the selected search type. Please ensure that the input directory contains only .fa files (.fai files are accepted and automatically skipped from the search)')
        sys.exit()

    # run searchBruteForce
    print("Search START")
    start_time = time.time()
    subprocess.run([corrected_origin_path+"searchBruteForce", str(genomeDir), str(
        filePAM),	str(fileGuide), str(mm), str(result), str(th), str(r), str(variant)])
    print("Search END")
    print("Search runtime: %s seconds" % (time.time() - start_time))

    if "-scores" in sys.argv[1:]:
        try:
            idx_genome_fasta = (sys.argv).index("-scores") + 1
            idx_genome_fasta = os.path.realpath(sys.argv[idx_genome_fasta])
        except:
            print(
                "ERROR! Please select the directory containing the fasta files of the genome")
            sys.exit()

        pam_guide = len(open(filePAM).readline().split(" ")[0])
        pam_at_beginning = int(open(filePAM).readline().split(" ")[1])
        if (pam_guide != 23 or pam_at_beginning < 0):  # Also block scoring pam at beginning
            print("WARNING: The CFD score and the Doench score can be calculated only for guides with 20bp and a 3bp PAM (SpCas9)")

        target_filename = os.path.realpath(result)
        subprocess.run([corrected_origin_path+'Python_Scripts/Scores/scores.py', target_filename +
                        '.targets.txt', idx_genome_fasta + "/", str(filePAM), str(fileGuide)])


def scores():
    if (len(sys.argv) < 6 or 'help' in sys.argv[1:]):
        print("WARNING: Too few arguments to function annotate-results. Please provide:\n",
              "\nEXAMPLE CALL: crispritz.py scores resultsFile.txt pamFile.txt guideFiles.txt genomeDirectoryInFastaFormat\n",
              "\n<resultsFile>: Targets file containing all genomic targets for the guides set",
              "\n<pamFile>: Text file containing the PAM sequence (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence",
              "\n<guidesFile>: Text file containing one or more guides (including a number of Ns equal to the length of the PAM sequence)",
              "\n<genomeDirectoryInFastaFormat>: Directory containing the genome in .fa or .fasta format, necessary to extract sequences for Doench Score Function",
              )
        sys.exit()

    resultFile = os.path.realpath(sys.argv[2])
    filePAM = os.path.realpath(sys.argv[3])
    fileGuide = os.path.realpath(sys.argv[4])
    genomeDir = os.path.realpath(sys.argv[5])+"/"

    checkExistance(resultFile, 'f')
    checkExistance(filePAM, 'f')
    checkExistance(fileGuide, 'f')
    checkExistance(genomeDir, 'd')

    pam_guide = len(open(filePAM).readline().split(" ")[0])
    pam_at_beginning = int(open(filePAM).readline().split(" ")[1])

    if (pam_guide != 23 or pam_at_beginning < 0):  # Also block scoring pam at beginning
        print("WARNING: The CFD score and the Doench score can be calculated only for guides with 20bp and a 3bp PAM (SpCas9)")
        # sys.exit()
    subprocess.run([corrected_origin_path+'Python_Scripts/Scores/scores.py',
                    resultFile, genomeDir, str(filePAM), str(fileGuide)])


def annotateResults():
    if (len(sys.argv) < 5 or 'help' in sys.argv[1:]):  # was 6
        print("WARNING: Too few arguments to function annotate-results. Please provide:\n",
              "\nEXAMPLE CALL: crispritz.py annotate-results resultsFile.txt annotationsFile.bed outputFile\n",
              "\n<resultsFile>: Targets file containing all genomic targets for the guides set",
              # Bed file containing annotation
              "\n<annotationsFile>: Text file containing the annotations in .bed format",
              "\n<outputFile>: Name of output file",
              "\n--change-ID <sampleIDfile> : (Optional) Change the samples, population and superpopulation IDs. DEFAULT: the default IDs are taken from the 1000 genome project (used for Human Genome hg19 and hg38)"
              )
        sys.exit()

    #guidesFile = os.path.realpath(sys.argv[2])
    resultsFile = os.path.realpath(sys.argv[2])
    annotationsFile = os.path.realpath(sys.argv[3])
    outputFile = os.path.realpath(sys.argv[4])

    checkExistance(resultsFile, 'f')
    checkExistance(annotationsFile, 'f')

    sampleIDfile = corrected_origin_path + \
        'Python_Scripts/ProcessData/samples_1000genomeproject.txt'
    if "--change-ID" in sys.argv[1:]:
        sampleIDfile = (sys.argv).index("--change-ID") + 1
        try:
            sampleIDfile = os.path.realpath(sys.argv[sampleIDfile])
        except:
            print('Warning! Missing sampleID file')
            sys.exit()
    checkExistance(sampleIDfile, 'f')
# NOTE sys.argv[-1] can be either sampleIDfile or 'Step [x/y] - Annotation', and it's just needed for having an output of the state of the
    # analysis for the user
    step = sys.argv[-1]
    if 'Step' not in step:
        print("Annotation START")
    start_time = time.time()
    subprocess.run([corrected_origin_path+'Python_Scripts/Annotator/annotator.py',
                    annotationsFile, resultsFile, outputFile, sampleIDfile, sys.argv[-1]])
    if 'Step' not in step:
        print("Annotation END")
        print("Annotation runtime: %s seconds" % (time.time() - start_time))


def genomeEnrichment_subprocess_VCF(altfile, genfile, dirGenome, doit, dirVCFFiles):
    subprocess.run([corrected_origin_path + 'Python_Scripts/Enrichment/enricher.py',
                    altfile, genfile, dirGenome.split('/')[-1], str(doit), dirVCFFiles])


def genomeEnrichment():
    '''
    Enrich the genome by parsing the vcf files and then replacing the nucleotides in the reference genome sequence. Compatible with indels
    '''
    if (len(sys.argv) < 4 or 'help' in sys.argv[1:]):
        print("WARNING: Too few arguments to function add-variants. Please provide:",
              "\nEXAMPLE CALL: crispritz.py add-variants vcfFilesDirectory/ genomeDirectory/\n",
              "\n\n<vcfFilesDirectory> : Directory containing VCF files, need to be separated into single chromosome files (multi-sample files will be collapsed into one fake individual)",
              "\n\n<genomeDirectory> : Directory containing a genome in .fa or .fasta format, need to be separated into single chromosome files.")
        sys.exit()

    dirVCFFiles = os.path.realpath(sys.argv[2])
    dirGenome = os.path.realpath(sys.argv[3])
    doit = 'no'
    if len(sys.argv) > 4 and sys.argv[4] == 'true':
        doit = 'yes'
    checkExistance(dirVCFFiles, 'd')
    checkExistance(dirGenome, 'd')
    listChrs = os.listdir(dirVCFFiles)
    # listChrs = glob.glob(dirVCFFiles+'/*.vcf.gz')
    
    for file in listChrs:
        if file.endswith('.tbi'): #remove .tbi files
            listChrs.remove(file)
        # if 'tbi' in elem:  # remove .tbi files in vcf dir to avoid errors in exec
            # listChrs.remove(elem)

    # print(listChrs)

    chr_with_vcf = set()
    # VCF file must contain '.chrN.' with N = number or letter of chr
    for vcf_chr in [f for f in listdir(dirVCFFiles) if isfile(join(dirVCFFiles, f)) and 'vcf.gz' in f]:
        for x in vcf_chr.split('.'):
            if 'chr' in x:
                chr_with_vcf.add(x)
    list_file_ends = [f for f in listdir(
        dirGenome) if isfile(join(dirGenome, f))]
    file_ends = list_file_ends[0].split('.')[-1]  # For .fa or .fasta
    contains_enr = ''
    if '.enriched.' in list_file_ends[0]:
        contains_enr = '.enriched'
    elif '.indels.' in list_file_ends[0]:
        contains_enr = '.indels'
    chr_wihtout_vcf = set([f.split(contains_enr + '.fa')[0]
                           for f in listdir(dirGenome) if isfile(join(dirGenome, f))]) - chr_with_vcf

    # dirParsedFiles = "./parsed_vcf_files/"

    # if not (os.path.isdir(dirParsedFiles)):
    #     os.makedirs(dirParsedFiles)

    dirEnrichedGenome = "./variants_genome/"
    if not (os.path.isdir(dirEnrichedGenome)):
        os.makedirs(dirEnrichedGenome)

    os.chdir(dirEnrichedGenome)
    if not (os.path.isdir("./SNPs_genome/")):
        os.makedirs("./SNPs_genome/")
    if not (os.path.isdir("./INDELs_genome/")):
        os.makedirs("./INDELs_genome/")
    os.chdir("./INDELs_genome/")

    memo = open('change_version.txt', 'w')

    memo.write('CRISPRitz indels process is now obsolete and has been removed, if you want to process indels you can download our new tool CRISPRme, https://github.com/samuelecancellieri/CRISPRme'+'\n')
    memo.write('Thank you')

    memo.close()

    os.chdir('../')
    # os.chdir(dirParsedFiles)
    pool = multiprocessing.Pool(4)

    print("Variants Extraction and Processing START")
    start_time = time.time()
    for elem in listChrs:
        # os.chdir(dirParsedFiles) #enter in parsed file directory
        chrom = ""
        for cut in elem.split('.'):
            if 'chr' in cut:
                chrom = cut
        # subprocess.run([corrected_origin_path +
        #                'Python_Scripts/Enrichment/bcf_query.sh', dirVCFFiles+"/"+elem, chrom])
        #altfile = str(chrom + '.alt')
        altfile = dirVCFFiles+"/"+elem
        genfile = str(dirGenome+'/' + chrom + contains_enr + '.fa')
        # pool to process 4 vcf file in parallel
        pool.apply_async(genomeEnrichment_subprocess_VCF, args=(
            altfile, genfile, dirGenome, doit, dirVCFFiles))
    # wait until all threads are completed than join
    pool.close()
    pool.join()

    for f in chr_wihtout_vcf:  # Move chromosomes without vcf to enriched directory and change name adding '.enriched.'
        subprocess.run(['cp', dirGenome + '/' + f + contains_enr + '.' + file_ends, './SNPs_genome/' +
                        dirGenome.split('/')[-1] + '_enriched/' + f + '.enriched.' + file_ends])

    print("Variants Extraction and Processing END")
    print("Runtime: %s seconds" % (time.time() - start_time))

    # os.chdir("../")
    # listChrs = os.listdir(dirParsedFiles)
    # dirEnrichedGenome = "./variants_genome/"

    # if not (os.path.isdir(dirEnrichedGenome)):
    #     os.makedirs(dirEnrichedGenome)

    # os.chdir(dirEnrichedGenome)
    # if not (os.path.isdir("./SNPs_genome/")):
    #     os.makedirs("./SNPs_genome/")
    # if not (os.path.isdir("./INDELs_genome/")):
    #     os.makedirs("./INDELs_genome/")

    # print("Variants Adding START")

    # start_time = time.time()
    # for f in listChrs:
    #     splitf = f.split('.')
    #     altfile = str('../'+dirParsedFiles+splitf[0]+'.alt')
    #     genfile = str(dirGenome+'/'+splitf[0] + contains_enr + '.fa')
    #     if not isfile(genfile) and not isfile(str(dirGenome+'/'+splitf[0] + contains_enr + '.fasta')):
    #         continue
    #     print("Adding Variants to:", splitf[0])
    #     subprocess.run([corrected_origin_path + 'Python_Scripts/Enrichment/enricher.py',
    #                     altfile, genfile, dirGenome.split('/')[-1]])

    # for f in chr_wihtout_vcf:  # Move chromosomes without vcf to enriched directory and change name adding '.enriched.'
    #     subprocess.run(['cp', dirGenome + '/' + f + contains_enr + '.' + file_ends, './SNPs_genome/' +
    #                     dirGenome.split('/')[-1] + '_enriched/' + f + '.enriched.' + file_ends])
    #     # subprocess.run(['cp', dirGenome + '/' + f + contains_enr + '.' + file_ends,
    #     #                 './INDELs_genome/' + dirGenome.split('/')[-1] + '_enriched/' + f + '.indels.' + file_ends])

    # print("Variants Adding END")
    # print("Runtime: %s seconds" % (time.time() - start_time))

    # os.chdir("../")
    # shutil.rmtree(dirParsedFiles) DO NOT REMOVE ALT FILE TO REUSE THEM LATER


def generateReport():
    if (len(sys.argv) < 9 or 'help' in sys.argv[1:]):  # was 10
        print("WARNING: Too few arguments to function generate-report. Please provide:\n",
              "\nEXAMPLE CALL: crispritz.py GAGTCCGAGCAGAAGAAGAANNN -mm 4 -annotation annotationSummaryFile.txt -extprofile guideExtendedProfile.xls -gecko -sumref referenceAnnotationSummaryFile.txt\n",
              "\n<guide>: (Optional) A guide present in the analyzed set",
              "\n-mm <mm_num>: Number of mismatches to analyze",
              "\n-annotation <annotationSummaryFile>: Count files for genomic annotations",
              "\n-extprofile <guideExtendedProfile>: Extended profile file",
              "\n-gecko: (Optional) Tag to activate gecko dataset comparison",
              "\n-sumref <referenceAnnotationSummaryFile>: (Optional) Create a barplot comparing reference genome results with enriched genome results. If the <guide> option is used, the barplot will take into account only the targets found with that specific guide\n")
        sys.exit()

    # NOTE barplot is ok only when -annotator and profile, ext are from var, -sumref is from ref
    if '-mm' in sys.argv[2]:
        guide = 'no'
    else:
        guide = sys.argv[2]

    mm = 0
    if "-mm" in sys.argv[1:]:
        try:
            mm = (sys.argv).index("-mm") + 1
            mm = (sys.argv[mm])
            int(mm)
        except:
            print(
                "ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is a number)")
            sys.exit()

    extProfileFile = "no"
    if "-extprofile" in sys.argv[1:]:
        extProfileFile = (sys.argv).index("-extprofile") + 1
        extProfileFile = os.path.realpath(sys.argv[extProfileFile])
        checkExistance(extProfileFile, 'f')

    countdir = 'no'
    if '-annotation' in sys.argv[1:]:
        countdir = (sys.argv).index("-annotation") + 1
        countdir = os.path.realpath(sys.argv[countdir])
        checkExistance(countdir, 'f')

    summaryOne = "no"
    if "-sumref" in sys.argv[1:]:
        summaryOne = (sys.argv).index("-sumref") + 1
        summaryOne = os.path.realpath(sys.argv[summaryOne])
        checkExistance(summaryOne, 'f')

    summaryTwo = countdir

    geckoProfile = "no"

    if "-gecko" in sys.argv[1:]:
        geckoProfile = corrected_origin_path + \
            'Python_Scripts/Plot/gecko.summary.total.Annotation.summary.txt'

    if '-ws' in sys.argv[:]:
        web_server = '-ws'
    else:
        web_server = ''
    if '-sample' in sys.argv[:]:
        sample_name = sys.argv.index('-sample')
        sample_name = sys.argv[sample_name + 1]
        sample_opt = '-sample'
    else:
        sample_name = ''
        sample_opt = ''

    subprocess.run([corrected_origin_path + 'Python_Scripts/Plot/radar_chart.py', str(guide), str(mm), str(summaryTwo), str(extProfileFile), str(summaryOne),
                    str(geckoProfile), web_server, sample_opt, sample_name])


def removeFile(to_remove):
    '''
    Remove file given his path in input
    '''
    try:
        os.remove(to_remove)
    except:
        pass

# return version


def version():
    print("CRISPRitz v2.5.8")


# HELP FUNCTION
def callHelp():

    print("help:\n",
          "\nALL FASTA FILEs USED BY THE SOFTWARE MUST BE UNZIPPED AND CHROMOSOME SEPARATED, ALL VCFs USED BY THE SOFTWARE MUST BE ZIPPED AND CHROMOSOME SEPARATED",
          "\n"
          "\ncrispritz.py add-variants FUNCTION TO ADD VARIANTS DATA TO A FASTA GENOME",
          "\ncrispritz.py index-genome FUNCTION TO CREATE GENOME INDEX TO PERFORM FAST SEARCHES WITH BULGES",
          "\ncrispritz.py search FUNCTION TO PERFORM SEARCHES ON A GENOME (INDEXED OR PLAIN FASTA)",
          "\ncrispritz.py scores FUNCTION TO CALCULATE THE CFD SCORE FOR A LIST OF TARGETS",
          "\ncrispritz.py annotate-results FUNCTION TO ADD GENOMIC INFORMATION TO TARGETS RESULTS",
          "\ncrispritz.py generate-report FUNCTION TO GENERATE GRAPHICAL REPORT FOR A SPECIFIC GUIDE",
          "\n",
          "\nADD help TO ANY FUNCTION TO VISUALIZE A BRIEF HELP PAGE (example: crispritz.py index-genome help)\n")


if len(sys.argv) < 2:
    callHelp()
elif sys.argv[1] == "version":
    version()
elif sys.argv[1] == "index-genome":
    indexGenome()
elif sys.argv[1] == "search" and ("-bDNA" in sys.argv[1:] or "-bRNA" in sys.argv[1:]):
    searchTST()
elif sys.argv[1] == "search":
    searchBruteForce()
elif sys.argv[1] == "scores":
    scores()
elif sys.argv[1] == "add-variants":
    genomeEnrichment()
elif sys.argv[1] == "annotate-results":
    annotateResults()
elif sys.argv[1] == "generate-report":
    generateReport()
else:
    print("ERROR! \"" + sys.argv[1] + "\" is not an allowed!")

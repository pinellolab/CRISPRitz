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
# path where this file is located
origin_path = os.path.dirname(os.path.realpath(__file__))
#conda path
conda_path="opt/crispritz/"
#path corrected to use with conda
corrected_origin_path = origin_path[:-3]+conda_path

def checkExistance(f_path, element):
	if element == 'f':		#check file
		if not isfile(f_path):
			print('ERROR! ' + f_path + ' is not a file.')
			sys.exit()
	else:	#check directory
		if not isdir(f_path):
			print('ERROR! ' + f_path + ' is not a directory.')
			sys.exit()

#ALL ACTIVE FUNCTIONS IN CRISPRITZ
def indexGenome():
	if (len(sys.argv) < 7 or 'help' in sys.argv[1:]):
		print("WARNING: Too few arguments to function index-genome. Please provide:\n",
		"\nEXAMPLE CALL: crispritz.py index-genome name_genome genomeDirectory(FASTA)/ pamFile.txt -bMax 2\n",
		"\n<name_genome>: Name of the genome to create",
		"\n<genomeDirectory>: Directory containing a genome in .fa or .fasta format, need to be separated into single chromosome files.",
		"\n<pamFile>: Text file containing the PAM (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence",
		"\n-bMax <maxBulges_num>: Number of bulges allowed for the search phase",
		"\n-th <num_thread>: (Optional) Number of threads to use. Default uses half of the available threads")
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

	TSTgenome = PAM + "_" + max_bulges + "_" + nameGenome					# name of the genome in TST format
	dirTSTgenome = "./genome_library/" + TSTgenome		# dir of the genome in TST format

	print(TSTgenome, "Indexing generation:")

	# # variant
	# variant = 0
	# if "-var" in sys.argv[1:]:
	# 	variant = 1
	
	# read number of threads
	th = 0
	if "-th" in sys.argv[1:]:
		try:
			th = (sys.argv).index("-th") + 1
			th = int(sys.argv[th])
		except:
			print("ATTENTION! Check the thread option: -th <th_num> (th_num is an integer)")
			sys.exit()

	if os.path.isdir(dirTSTgenome):							# check if TSTgenome dir exists
		shutil.rmtree(dirTSTgenome)							# remove old TSTgenome dir
	os.makedirs(dirTSTgenome)								# build new TSTgenome dir
	os.chdir(dirTSTgenome)									# move into the TSTgenome dir

	# run buildTST
	start_time = time.time()
	for f in listChrs:
		print("Indexing:", f)
		subprocess.run([corrected_origin_path+"buildTST",str(dirGenome)+"/"+str(f), str(dirPAM), str(th), max_bulges])
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
			print("ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is a number)")
			sys.exit()

	# read number of bulge RNA
	bRNA = 0
	if "-bRNA" in sys.argv[1:]:
		try:
			bRNA = (sys.argv).index("-bRNA") + 1
			bRNA = int(sys.argv[bRNA])
		except:
			print("ATTENTION! Check the bulge RNA option: -bRNA <bRNA_num> (bRNA_num is a number)")
			sys.exit()

	# read number of bulge DNA
	bDNA = 0
	if "-bDNA" in sys.argv[1:]:
		try:
			bDNA = (sys.argv).index("-bDNA") + 1
			bDNA = int(sys.argv[bDNA])
		except:
			print("ATTENTION! Check the bulge DNA option: -bDNA <bDNA_num> (bDNA_num is a number)")
			sys.exit()

	if (int(max_bulges) < bDNA):
		print("WARNING! Max available bulges (" + max_bulges + ") is smaller than input DNA bulges (" + str(bDNA) + "). Using " + max_bulges + " bulges")
		bDNA = max_bulges
	if (int(max_bulges) < bRNA):
		print("WARNING! Max available bulges (" + max_bulges + ") is smaller than input RNA bulges (" + str(bRNA) + "). Using " + max_bulges + " bulges")
		bRNA = max_bulges
	# read number of threads
	th = 1
	if "-th" in sys.argv[1:]:
		try:
			th = (sys.argv).index("-th") + 1
			th = int(sys.argv[th])
		except:
			print("ATTENTION! Check the thread option: -th <th_num> (th_num is an integer)")
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
			print("ERROR! Please select the directory containing the fasta files of the genome")
			sys.exit()
		checkExistance(idx_genome_fasta, 'd')
	# Check input correctness
	file_correct_ext = [f for f in listdir(dirTSTgenome) if isfile(join(dirTSTgenome, f)) and not f.endswith('.bin')]	#Get files not ending with .bin
	if len(file_correct_ext) != 0:		# Some files do not have .bin
		print('ERROR! The directory contains files that are not compatible with the selected search type. Please ensure that the input directory contains only .bin files')
		sys.exit()
	
	# run searchOnTST
	print("Search START")
	start_time = time.time()
	subprocess.run([corrected_origin_path+"searchTST", str(dirTSTgenome), str(fileGuide), str(mm), str(bDNA), str(bRNA), str(PAM), str(nameResult), str(r), str(th), max_bulges])
	print("Search END")
	print("Search runtime: %s seconds" % (time.time() - start_time))
	if "-scores" in sys.argv[1:]:

		try:
			idx_genome_fasta = (sys.argv).index("-scores") + 1
			idx_genome_fasta = os.path.realpath(sys.argv[idx_genome_fasta])
		except:
			print("ERROR! Please select the directory containing the fasta files of the genome")
			sys.exit()
		
		# pam_len = int(open(PAM).readline().split(" ")[1])
		# if (pam_len < 0):
		# 	pam_begin = True
		# else:
		# 	pam_begin = False
		pam_guide = len(open(PAM).readline().split(" ")[0])
		pam_at_beginning = int(open(PAM).readline().split(" ")[1])
		if (pam_guide != 23 or pam_at_beginning < 0): #Also block scoring pam at beginning
			print("WARNING: The CFD score and the Doench score can be calculated only for guides with 20bp and a 3bp PAM")
			sys.exit()
		# pam_seq_check_ngg = open(PAM).readline().split(" ")[0].upper()
		# if ("NGG" not in pam_seq_check_ngg):
		# 	# print("WARNING: The model used for the CFD and Doench scores are based on the NGG PAM, the scores may not be valid for other PAMs")
		# 	print("WARNING: The model used for the CFD and Doench scores are based on the NGG PAM")
		# 	sys.exit()
		target_filename = os.path.realpath(nameResult)
		subprocess.run([corrected_origin_path+'Python_Scripts/Scores/scores.py', target_filename + '.targets.txt', idx_genome_fasta + "/", str(PAM), str(fileGuide)])

def searchBruteForce():
	if (len(sys.argv) < 9 or 'help' in sys.argv[1:]):
		print("WARNING: Too few arguments to function search. Please provide:\n",
		"\nEXAMPLE CALL MISMATCHES ONLY: crispritz.py search genomeDirectory/ pamFile.txt guidesFile.txt outputFile -mm 4 -th 4 -var -scores genomeDirectory(FASTA)/ -t\n",
		"\nEXAMPLE CALL MISMATCHES + BULGES: crispritz.py search indexGenomeDirectory/ pamFile.txt guidesFile.txt outputFile -mm 4 -bDNA 1 -bRNA 2 -th 4 -var -scores genomeDirectory(FASTA)/ -t\n",
		"\n<genomeDirectory>: Directory containing a genome in .fa or .fasta format (.bin format if bulges present), need to be separated into single chromosome files",
		"\n<pamFile>: Text file containing the PAM sequence (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence" ,
		"\n<guidesFile>: Text file containing one or more guides (including a number of Ns equal to the length of the PAM sequence)", 
		"\n<outputFile>: Name of output file",
		"\n-mm <mm_num>: Number of allowed mismatches",
		"\n-bRNA <bRNA_num>: (Optional) Size of RNA bulges",
		"\n-bDNA <bDNA_num>: (Optional) Size of DNA bulges" ,
		"\n-th <num_thread>: (Optional) Number of threads to use",
		"\n-var: (Optional) Tag to activate search with IUPAC nomenclature",
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
			print("ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is an integer)")
			sys.exit()

	# read number of mismatches
	th = 10000000
	if "-th" in sys.argv[1:]:
		try:
			th = (sys.argv).index("-th") + 1
			th = int(sys.argv[th])
		except:
			print("ATTENTION! Check the thread option: -th <th_num> (th_num is an integer)")
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
	variant = 0
	if "-var" in sys.argv[1:]:
		variant = 1
	
	# Check '-scores' directory
	if "-scores" in sys.argv[1:]:
		try:
			idx_genome_fasta = (sys.argv).index("-scores") + 1
			idx_genome_fasta = os.path.realpath(sys.argv[idx_genome_fasta])
		except:
			print("ERROR! Please select the directory containing the fasta files of the genome")
			sys.exit()
		checkExistance(idx_genome_fasta, 'd')

	# Check input correctness
	checkExistance(genomeDir, 'd')
	checkExistance(filePAM, 'f')
	checkExistance(fileGuide, 'f')
	
	file_correct_ext = [f for f in listdir(genomeDir) if isfile(join(genomeDir, f)) and not (f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.fai'))]
	
	if len(file_correct_ext) != 0:		#Some file other than .fa, .fasta or .fai are present in the directory
		print('ERROR! The directory contains files that are not compatible with the selected search type. Please ensure that the input directory contains only .fa files (.fai files are accepted and automatically skipped from the search)')
		sys.exit()

	# run searchBruteForce
	print("Search START")
	start_time = time.time()
	subprocess.run([corrected_origin_path+"searchBruteForce", str(genomeDir), str(filePAM),	str(fileGuide), str(mm), str(result), str(th), str(r), str(variant)])
	print("Search END")
	print("Search runtime: %s seconds" % (time.time() - start_time))

	if "-scores" in sys.argv[1:]:
		try:
			idx_genome_fasta = (sys.argv).index("-scores") + 1
			idx_genome_fasta = os.path.realpath(sys.argv[idx_genome_fasta])
		except:
			print("ERROR! Please select the directory containing the fasta files of the genome")
			sys.exit()
		
		# pam_len = int(open(filePAM).readline().split(" ")[1])
		# if (pam_len < 0):
		# 	pam_begin = True
		# else:
		# 	pam_begin = False
		 
		pam_guide = len(open(filePAM).readline().split(" ")[0])
		pam_at_beginning = int(open(filePAM).readline().split(" ")[1])
		if (pam_guide != 23 or pam_at_beginning < 0): #Also block scoring pam at beginning
			print("WARNING: The CFD score and the Doench score can be calculated only for guides with 20bp and a 3bp PAM")
			sys.exit()
		# pam_seq_check_ngg = open(filePAM).readline().split(" ")[0].upper()
		# if ("NGG" not in pam_seq_check_ngg):
		# 	# print("WARNING: The model used for the CFD and Doench scores are based on the NGG PAM, the scores may not be valid for other PAMs")
		# 	print("WARNING: The model used for the CFD and Doench scores are based on the NGG PAM")
		# 	sys.exit()
		target_filename = os.path.realpath(result)
		subprocess.run([corrected_origin_path+'Python_Scripts/Scores/scores.py', target_filename + '.targets.txt', idx_genome_fasta + "/", str(filePAM), str(fileGuide)])

def scores():
	if (len(sys.argv) < 6 or 'help' in sys.argv[1:]):
		print("WARNING: Too few arguments to function annotate-results. Please provide:\n",
			"\nEXAMPLE CALL: crispritz.py scores resultsFile.txt genomeDirectoryInFastaFormat\n",
			"\n<resultsFile>: Targets file containing all genomic targets for the guides set",
			"\n<pamFile>: Text file containing the PAM sequence (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence" ,
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

	if (pam_guide != 23 or pam_at_beginning < 0):	#Also block scoring pam at beginning
		print("WARNING: The CFD score and the Doench score can be calculated only for guides with 20bp and a 3bp PAM")
		sys.exit()
	subprocess.run([corrected_origin_path+'Python_Scripts/Scores/scores.py', resultFile, genomeDir, str(filePAM), str(fileGuide)])

def annotateResults():
	if (len(sys.argv) < 5 or 'help' in sys.argv[1:]):	#was 6
		print("WARNING: Too few arguments to function annotate-results. Please provide:\n",
		"\nEXAMPLE CALL: crispritz.py annotate-results resultsFile.txt annotationsFile.bed outputFile\n",
		#"\n<profileFile>: Profile File containing all information grouped by guide",
		"\n<resultsFile>: Targets file containing all genomic targets for the guides set",
		"\n<annotationsFile>: Text file containing the annotations in .bed format", # -> Bed file containing annotation
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

	sampleIDfile = corrected_origin_path + 'Python_Scripts/ProcessData/samples_1000genomeproject.txt'
	if "--change-ID" in sys.argv[1:]:
		sampleIDfile = (sys.argv).index("--change-ID") + 1
		try:
			sampleIDfile = os.path.realpath(sys.argv[sampleIDfile])
		except:
			print('Warning! Missing sampleID file')
			sys.exit()
	checkExistance(sampleIDfile, 'f')
    #NOTE sys.argv[-1] can be either sampleIDfile or 'Step [x/y] - Annotation', and it's just needed for having an output of the state of the 
	#analysis for the user
	step = sys.argv[-1]
	if 'Step' not in step:
		print("Annotation START")
	start_time = time.time()
	subprocess.run([corrected_origin_path+'Python_Scripts/Annotator/annotator.py', annotationsFile, resultsFile, outputFile, sampleIDfile ,sys.argv[-1]])
	if 'Step' not in step:
		print("Annotation END")
		print("Annotation runtime: %s seconds" % (time.time() - start_time))

def genomeEnrichmentWithBedtools():
	'''
	Enrich the genome using the bedtools program. Not used since genomeEnrichment() is now compatible with indels
	'''
	if (len(sys.argv) < 4 or 'help' in sys.argv[1:]):
		print("WARNING: Too few arguments to function add-variants. Please provide:\n",
		"\nEXAMPLE CALL: crispritz.py add-variants vcfFilesDirectory/ genomeDirectory/\n",
		"\n<vcfFilesDirectory>: Directory containing VCF files, MUST BE separated into single chromosome files",
		"\n<genomeDirectory>: Directory containing a genome in .fa or .fasta format, MUST BE separated into single chromosome files\n")
		sys.exit()
	dirVCFFiles = os.path.realpath(sys.argv[2])		#VCF files must contain .chrN. where N is the chr name in their filename
	dirGenome = os.path.realpath(sys.argv[3])		#Chromosomes filename must be chrN
	listChrs = os.listdir(dirVCFFiles)

	checkExistance(dirVCFFiles, 'd')
	checkExistance(dirGenome, 'd')
	genomeName=dirGenome.split('/')
	genomeName=genomeName[-1]

	vcfList = glob.glob(dirVCFFiles+'/*.vcf.gz')

	print('Variant Extraction START')
	start_time = time.time()

	cmds_list = [[corrected_origin_path+'Python_Scripts/Enrichment/index_vcf.sh', file_name] for file_name in vcfList]
	procs_list = [Popen(cmd, stdout=PIPE, stderr=PIPE) for cmd in cmds_list]

	for proc in procs_list:
		proc.wait()

	print('Variant Extraction END')
	print("Runtime: %s seconds" % (time.time() - start_time))

	dirEnrichedGenome = genomeName+"_variants"

	if (os.path.isdir(dirEnrichedGenome)):
		shutil.rmtree(dirEnrichedGenome)
		os.makedirs(dirEnrichedGenome)
	else:
		os.makedirs(dirEnrichedGenome)
	
	os.chdir(dirEnrichedGenome)

	chrList = []
	vcfNameList = []

	print("Variants Adding START")
	start_time = time.time()

	for vcf_file in vcfList:
		if vcf_file.endswith('.csi'):
			continue
		vcfName = vcf_file.split('/')[-1]
		vcfNameList.append(vcfName)
		chr_name_start = vcf_file.find('.chr')
		chr_name_end = vcf_file.find('.', chr_name_start + 1)
		chr_name = vcf_file[chr_name_start + 1 : chr_name_end]
		chrList.append(chr_name)

	cmds_list = [[corrected_origin_path+'Python_Scripts/Enrichment/add_variants.sh', dirVCFFiles+'/'+vcfNameList[count], dirGenome + '/' + chrName + '.fa', chrName] for count,chrName in enumerate(chrList)]
	procs_list = [Popen(cmd, stdout=PIPE, stderr=PIPE) for cmd in cmds_list]

	for proc in procs_list:
		proc.wait()

	print("Variants Adding END")
	print("Runtime: %s seconds" % (time.time() - start_time))

	os.chdir("../..")

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
	checkExistance(dirVCFFiles, 'd')
	checkExistance(dirGenome, 'd')
	listChrs = os.listdir(dirVCFFiles)
	chr_with_vcf = set()
	for vcf_chr in [f for f in listdir(dirVCFFiles) if isfile(join(dirVCFFiles, f))]:	#VCF file must contain '.chrN.' with N = number or letter of chr
		for x in vcf_chr.split('.'):
			if 'chr' in x:
				chr_with_vcf.add(x)
	list_file_ends = [f for f in listdir(dirGenome) if isfile(join(dirGenome, f))]
	file_ends = list_file_ends[0].split('.')[-1]	#For .fa or .fasta
	contains_enr = ''
	if '.enriched.' in list_file_ends[0]:
		contains_enr = '.enriched'
	elif '.indels.' in list_file_ends[0]:
		contains_enr = '.indels'
	chr_wihtout_vcf = set([f.split(contains_enr + '.fa')[0] for f in listdir(dirGenome) if isfile(join(dirGenome, f))]) - chr_with_vcf

	dirParsedFiles = "./parsed_vcf_files/"

	if not (os.path.isdir(dirParsedFiles)):
		os.makedirs(dirParsedFiles)

	os.chdir(dirParsedFiles)

	print("Variants Extraction START")
	start_time = time.time()
	subprocess.run([corrected_origin_path + 'Python_Scripts/Enrichment/bcf_query.sh', dirVCFFiles+'/']) 
	print("Variants Extraction END")
	print("Runtime: %s seconds" % (time.time() - start_time))


	os.chdir("../")
	listChrs = os.listdir(dirParsedFiles)
	dirEnrichedGenome = "./variants_genome/"

	if not (os.path.isdir(dirEnrichedGenome)):
		os.makedirs(dirEnrichedGenome)

	os.chdir(dirEnrichedGenome)
	if not (os.path.isdir("./SNPs_genome/")):
		os.makedirs("./SNPs_genome/")
	if not (os.path.isdir("./INDELs_genome/")):
		os.makedirs("./INDELs_genome/")

	print("Variants Adding START")
	
	start_time = time.time()
	for f in listChrs:
		splitf = f.split('.')
		altfile = str('../'+dirParsedFiles+splitf[0]+'.alt')
		genfile = str(dirGenome+'/'+splitf[0] + contains_enr +'.fa')
		if not isfile(genfile) and not isfile (str(dirGenome+'/'+splitf[0] + contains_enr +'.fasta')):
			continue
		print("Adding Variants to:", splitf[0])
		subprocess.run([corrected_origin_path + 'Python_Scripts/Enrichment/enricher.py', altfile, genfile, dirGenome.split('/')[-1]]) 
	
	for f in chr_wihtout_vcf:		#Move chromosomes without vcf to enriched directory and change name adding '.enriched.'
		subprocess.run(['cp', dirGenome + '/' + f + contains_enr + '.' + file_ends, './SNPs_genome/' + dirGenome.split('/')[-1] + '_enriched/' + f + '.enriched.' + file_ends])
		subprocess.run(['cp', dirGenome + '/' + f + contains_enr + '.' + file_ends, './INDELs_genome/' + dirGenome.split('/')[-1] + '_enriched/' + f + '.indels.' + file_ends])

	print("Variants Adding END")
	print("Runtime: %s seconds" % (time.time() - start_time))

	os.chdir("../")
	shutil.rmtree(dirParsedFiles)

def generateReport():
	if (len(sys.argv) < 9 or 'help' in sys.argv[1:]):	#was 10
		print("WARNING: Too few arguments to function generate-report. Please provide:\n",
		"\nEXAMPLE CALL: crispritz.py GAGTCCGAGCAGAAGAAGAANNN -mm 4 -annotation annotationSummaryFile.txt -extprofile guideExtendedProfile.xls -gecko -sumref referenceAnnotationSummaryFile.txt\n",
		"\n<guide>: (Optional) A guide present in the analyzed set",
		"\n-mm <mm_num>: Number of mismatches to analyze",
		"\n-annotation <annotationSummaryFile>: Count files for genomic annotations",
		"\n-extprofile <guideExtendedProfile>: Extended profile file",
		#"\n-annotation <annotationFile>: Count files for genomic annotations",
		"\n-gecko: (Optional) Tag to activate gecko dataset comparison",
		"\n-sumref <referenceAnnotationSummaryFile>: (Optional) Create a barplot comparing reference genome results with enriched genome results. If the <guide> option is used, the barplot will take into account only the targets found with that specific guide\n")
		sys.exit()
		
	
	#NOTE barplot is ok only when -annotator and profile, ext are from var, -sumref is from ref
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
			print("ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is a number)")
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
	# geckoExonsCount = "no"
	# geckoIntronsCount = "no"
	# geckoPromotersCount = "no"
	# geckoDNAseCount = "no"
	# geckoCTCFCount = "no"

	if "-gecko" in sys.argv[1:]:
		geckoProfile = corrected_origin_path + \
				'Python_Scripts/Plot/gecko.summary.total.Annotation.summary.txt'
		# geckoProfile = corrected_origin_path + \
		# 		'Python_Scripts/Plot/gecko/gecko.reference.profile.xls'
		# geckoExonsCount = corrected_origin_path + \
		# 		'Python_Scripts/Plot/gecko/gecko.Exons.Count.txt'
		# geckoIntronsCount = corrected_origin_path + \
		# 		'Python_Scripts/Plot/gecko/gecko.Introns.Count.txt'
		# geckoPromotersCount = corrected_origin_path + \
		# 		'Python_Scripts/Plot/gecko/gecko.Promoters.Count.txt'
		# geckoDNAseCount = corrected_origin_path + \
		# 		'Python_Scripts/Plot/gecko/gecko.DNAse.Count.txt'
		# geckoCTCFCount = corrected_origin_path + \
		# 		'Python_Scripts/Plot/gecko/gecko.CTCF.Count.txt'
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
	# subprocess.run([corrected_origin_path + 'Python_Scripts/Plot/radar_chart.py', str(guide), str(mm), str(profileFile), str(extProfileFile), str(countdir), str(summaryOne), str(summaryTwo),
	# 				str(geckoProfile), str(geckoExonsCount), str(geckoIntronsCount), str(geckoPromotersCount), str(geckoDNAseCount), str(geckoCTCFCount)])
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

def processData():
	'''
	Function used to post-process the obtained data. Input files are the targets file obtained from the searches done one the reference and variant genomes.
	For the REF search, the targets are clustered and the top_1 are annotated, producing summaries for guide and positions.
	For VAR and VAR/REF search, the targets are also associated with events of pam creation, and the corresponding samples are written for each
	target.
	The -sample option is needed to indicate the directory containing the dictionaries for the sample extraction 
	The option --sample-create creates the dictionaries from the vcf directory given in input (heavy task) and then executes the -sample option.
	Example call:
	-For REF
	crispritz.py process-data -reftarget REFsearch.targets.txt pamfile.txt guidesfile.txt hg38Annotations.bed output.name -refgenome hg38_ref
	-For VAR
	crispritz.py process-data -vartarget VARsearch.targets.txt pamfile.txt guidesfile.txt hg38Annotations.bed output.name -sample dictionary_directory -refgenome hg38_var
	-For VAR/REF
	crispritz.py process-data -reftarget REFsearch.targets.txt -vartarget VARsearch.targets.txt pamfile.txt guidesfile.txt hg38Annotations.bed output.name -sample dictionary_directory -refgenome hg38_ref
	To create dictionary, change -sample dictionary_directory to --sample-create vcf_directory
	'''
	if (len(sys.argv) < 9 or 'help' in sys.argv[1:]):
		print("WARNING: Too few arguments to function process-data. Please provide:",
		"\nEXAMPLE CALL: crispritz.py process-data -reftarget resultsReferenceTargets.txt -vartarget resultsVariantTargets.txt pamFile.txt guidesFile.txt annotationsFile.bed outputFile -sample dictionaryDirectory/ -refgenome genomeReferenceDirectory/\n",
		"\n-reftarget <refTargetsFile> : Targets file, containing all genomic targets found in the reference genome",
		"\n-vartarget <varTargetsFile> : Targets file, containing all genomic targets found in the variant genome",
		"\n<pamFile> : Text file containing the PAM sequence (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence",
		"\n<guidesFile>: Text file containing one or more guides (including a number of Ns equal to the length of the PAM sequence)",
		"\n<annotationsFile> : Text file containing the annotations in .bed format",
		"\n<outputFile> : Name of output file",
		"\n-sample <dictionariesDirectory> : Add samples information and annotations to the generated file. Must be provided if -vartarget option is selected",
		"\n--sample-create <vcfFileDirectory> : (Optional) Creates samples dictionaries starting from the variant files (Heavy task), and add samples information to the generated file.",
	    "\n-refgenome <genomeDirectory> : Directory of reference Genome",
		"\n--change-ID <sampleIDfile> : (Optional) Change the samples, population and superpopulation IDs. DEFAULT: the default IDs are taken from the 1000 genome project (used for Human Genome hg19 and hg38)"
		)
		sys.exit()

	last_target_pos = 0 #Since pam is after -ref or -var, get the last one of them to know the pamfile position
	ref_result = "no"
	if "-reftarget" in sys.argv[1:]:
		ref_result = (sys.argv).index("-reftarget") + 1
		if ref_result > last_target_pos:
			last_target_pos = ref_result
		ref_result = os.path.realpath(sys.argv[ref_result])
		checkExistance(ref_result, 'f')
	
	var_result = "no"
	if "-vartarget" in sys.argv[1:]:
		var_result = (sys.argv).index("-vartarget") + 1
		if var_result > last_target_pos:
			last_target_pos = var_result
		var_result = os.path.realpath(sys.argv[var_result])
		checkExistance(var_result, 'f')
	if last_target_pos == 0:
		print("Warning! Please provide al least one target file (-reftarget or -vartarget options)")
		sys.exit()

	pam = 'no'
	try:
		pam = os.path.realpath(sys.argv[last_target_pos + 1])
	except:
		print('Warning! Missing pam file')
		sys.exit()
	checkExistance(pam, 'f')
	
	guides = 'no'
	try:
		guides = os.path.realpath(sys.argv[last_target_pos + 2])
	except:
		print('Warning! Missing guides file')
		sys.exit()
	checkExistance(guides, 'f')
	annotationFile = 'no'
	try:
		annotationFile = os.path.realpath(sys.argv[last_target_pos + 3])
	except:
		print('Warning! Missing annotation file')
		sys.exit()
	checkExistance(annotationFile, 'f')

	result = ''
	try:
		result = os.path.realpath(sys.argv[last_target_pos + 4])
	except:
		print('Warning! Missing output filename')
		sys.exit()
	
	sampleIDfile = corrected_origin_path + 'Python_Scripts/ProcessData/samples_1000genomeproject.txt'
	if "--change-ID" in sys.argv[1:]:
		sampleIDfile = (sys.argv).index("--change-ID") + 1
		try:
			sampleIDfile = os.path.realpath(sys.argv[sampleIDfile])
		except:
			print('Warning! Missing sampleID file')
			sys.exit()
		checkExistance(sampleIDfile, 'f')

	samples = False
	create_dict = False
	vcf_file_directory = 'no'
	if "--sample-create" in sys.argv[1:]:
		if "-vartarget" not in sys.argv[1:]:
			print('Warning! The sample analysis can only be performed when the -vartarget option is selected')
		vcf_file_directory = (sys.argv).index("--sample-create") + 1
		vcf_file_directory = os.path.realpath(sys.argv[vcf_file_directory])
		checkExistance(vcf_file_directory, 'd')
		samples = True
		create_dict = True
	
	dict_directory = 'no'
	if "-sample" in sys.argv[1:]:
		if "-vartarget" not in sys.argv[1:]:
			print('Warning! The sample analysis can only be performed when the -vartarget option is selected')
		dict_directory = (sys.argv).index("-sample") + 1
		dict_directory = os.path.realpath(sys.argv[dict_directory])
		checkExistance(dict_directory, 'd')
		samples = True
	
	if var_result != 'no' and ("--sample-create" not in sys.argv[1:] and  "-sample" not in sys.argv[1:] ): 	#-varresult selected but no dictionary provided
		print('Warning! If the option -vartarget is selected, the option -sample <dictionariesDirectory> must be provided')
		sys.exit()

	#Check for -refgenome
	reference_genome_dir = 'no'
	if "-refgenome" in sys.argv[1:]:
		reference_genome_dir = (sys.argv).index("-refgenome") + 1
		try:
			reference_genome_dir = os.path.realpath(sys.argv[reference_genome_dir])
		except:
			print('Warning! Missing <genomeDir> value for option -refgenome')
			sys.exit()
		checkExistance(reference_genome_dir, 'd')
	else:
		print('Error! Missing -refgenome option')
		sys.exit()
	
	process = subprocess.Popen(['wc', '-l', guides], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	num_guides = int(out.decode('UTF-8').split(' ')[0])
	
	#Create dictionaries
	if create_dict:
		print('Warning! Option --sample-create was selected. CRISPRitz will generate dictionaries that will be used to perform the sample analysis.\n', 
		'This heavy task will require about 20 minutes for each vcf file in the input directory, and will require about 80 Gb of space. Only SNP will be extracted from the vcf files.\n',
		'Please verify that the fasta files and the associated \"#CHROM\" column of the vcf file have the same name', sep = '')
		overwrite_dict = 'yes'
		dict_dir_name = vcf_file_directory.split('/')[-1]
		if not dict_dir_name:
			dict_dir_name = vcf_file_directory.split('/')[-2]
		if not os.path.isdir('dictionary_' + dict_dir_name):
			os.mkdir('dictionary_' + dict_dir_name)
		else:
			print('Directory', 'dictionary_' + dict_dir_name, 'already exists. Do you want to overwrite it? (yes/no)')
			overwrite_dict = input()
			while True:
				if overwrite_dict == 'yes':
					break
				elif overwrite_dict == 'no':
					break
				else:
					print('Please insert yes or no')
					overwrite_dict = input()
		if overwrite_dict == 'yes':
			vcf_files = [join(vcf_file_directory, f) for f in listdir(vcf_file_directory) if isfile(join(vcf_file_directory, f))]
			for vcf in vcf_files:
				dictionary_name = vcf[vcf.index('.chr') +1 : ]			# 'chr21.something'
				dictionary_name = dictionary_name[: dictionary_name.index('.')]
				print('Creating my_dict_' + dictionary_name + '.json')
				subprocess.run([corrected_origin_path+'Python_Scripts/ProcessData/dictionary_creation.py', vcf, dictionary_name, 'dictionary_' + dict_dir_name ])
			dict_directory = os.path.realpath('dictionary_' + dict_dir_name) #UPDATE dict_directory
		else:
			if '-vartarget' in sys.argv[:]:
				print('Warning! If the option -vartarget is selected, a dictionary must be provided.',
				'Please provide a dictionary directory using the option -sample <dictionariesDirectory> or create a new dictionary directory using the option --sample-create <vcfFileDirectory>', sep = '\n')
				exit()
	start_time_processdata = time.time()
	print('Processing Data START')
	if ref_result != 'no' and var_result == 'no':		#Process data for REFERENCE only result
		#Clustering
		print('Step [1/5] - Clustering Targets', end = '\r')
		step = 'Step [1/5] - Clustering Targets'
		start_time = time.time()
		subprocess.run([corrected_origin_path+'Python_Scripts/ProcessData/cluster.dict.py', ref_result, 'no', 'True', 'False', guides, result, step])
		print('Step [1/5] - Clustering Targets DONE: %.2f seconds' % round(time.time() - start_time,2))

		#Extract Top1
		print('Step [2/5] - Top1 Extraction', end = '\r')
		start_time = time.time()
		subprocess.run([corrected_origin_path+'Python_Scripts/ProcessData/extract_top.py', result + '.cluster.txt', result])
		#Get max mismatch, bulge dna and rna
		max_values = subprocess.Popen(['awk', 'BEGIN{m_m=0; b_d=0; b_r=0} NR>1 {if ($8>m_m) m_m=$8; if ($1=="DNA" && $9>b_d) b_d=$9; else if ($1=="RNA" && $9>b_r) b_r=$9} END{print m_m","b_d","b_r}', result + '.top_1.txt' ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = max_values.communicate()
		max_values = out.decode('UTF-8').strip().split(',')	#Get number of max mm, dna and rna bulge
		max_bulge = max_values[1]
		if int(max_values[2]) > int(max_values[1]):
			max_bulge = max_values[2]
		#Scores
		subprocess.call([corrected_origin_path+'Python_Scripts/Scores/scores.py', result + '.top_1.txt', reference_genome_dir, pam, guides])
		subprocess.call(['mv', result + '.top_1.txt.scores.txt', result + '.scores.txt'])
		subprocess.call(['mv', result + '.top_1.txt.targets.CFD.txt', result + '.targets.CFD.txt'])
		print('Step [2/5] - Top1 Extraction DONE: %.2f seconds' % round(time.time() - start_time,2))

		#Annotation of REF file
		print('Step [3/5] - Annotation', end = '\r')
		step = 'Step [3/5] - Annotation'
		start_time = time.time()
		subprocess.run(['crispritz.py','annotate-results', result + '.top_1.txt', annotationFile, result, step])
		print('\033[KStep [3/5] - Annotation DONE: %.2f seconds' % round(time.time() - start_time,2))

		#Summary by guide - pos 
		print('Step [4/5] - Creating Summaries : Summary by Guide and by Position', end = '\r')
		start_time = time.time()
		subprocess.run([corrected_origin_path+'Python_Scripts/ProcessData/summary_by_guide_position.py', result + '.cluster.txt', max_values[0], max_values[1], max_values[2], guides, result, 'No'])
		print('\033[KStep [4/5] - Creating Summaries DONE: %.2f seconds' % round(time.time() - start_time,2))

		#Create final file
		print('Step [5/5] - Generating Report', end = '\r')
		start_time = time.time()
		general_table = dict()	#GUIDE -> [target reference, off target reference]
		current_guide = 0
		with open(guides) as all_guides:
			for g in all_guides:
				g = g.strip()
				current_guide += 1
				print('\033Step [5/5] - Generating Report: Guide ' + str(current_guide) + '/' + str(num_guides), end = '\r')
				with open (result + '.targets.' + g + '.txt', 'w+') as final_result:
					subprocess.run(['head', '-1', result + '.Annotation.targets.txt'], stdout = final_result)	#Put header
				subprocess.run(['LC_ALL=C grep ' + g + ' ' + result + '.Annotation.targets.txt >> ' + result + '.targets.' + g + '.txt'], shell = True)	#NOTE shell = True
				# Generate file for general table
				general_table[g] = []
				df_profile = pd.read_csv(result + '.summary_by_guide.' + g + '.txt', sep = '\t')
				general_table[g].append(str(int(df_profile[(df_profile.Mismatches == 0) & (df_profile['Bulge Type'] == 'X')].iloc[0]['Targets in Reference'])))
				one_to_n_mms = []
				for i in range (1, int(max_values[0]) + 1 + int(max_bulge)):
					one_to_n_mms.append(sum(df_profile[((df_profile['Mismatches'] + df_profile['Bulge Size']) == i)]['Targets in Reference'].to_list()))
				general_table[g].append(str(int(sum(one_to_n_mms[1:]))) + ' (' + ' - '.join(str(int(x)) for x in one_to_n_mms) + ')')
		with open(result + '.general_target_count.txt', 'w+') as res_gen:
			res_gen.write('#Guide\tOn-Targets Reference\tOff-Targets Reference (' + ' - '.join([str(x) for x in range(int(max_values[0]) + 1 + int(max_bulge)) ]) + ' Mismatches + Bulges)\n')
			for g in general_table:
				res_gen.write(g + '\t' + '\t'.join(general_table[g]) + '\n')
		#Remove intermediate files
		removeFile(result + '.top_1.txt')
		removeFile(result + '.Annotation.targets.txt')
		print('Step [5/5] - Generating Report DONE: %.2f seconds' % round(time.time() - start_time,2))

	elif var_result != 'no':		#Process data for VARIANTS only result and VARIANTS/REFERENCE results
		file_to_cluster = result + '.total.txt'
		only_cluster = 'False'
		sum_guide_type = 'no'
		annotator_script_file = 'annotator_for_onlyvar.py'
		total_steps = '5'
		current_step = 1
		if ref_result != 'no':			#Steps for VARIANTS/REFERENCE analysis
			total_steps = '6'
			only_cluster = 'True'
			sum_guide_type = 'Uniq'
			annotator_script_file = 'annotator_cal_sample.py'
			#Extract Common, Semicommon and Unique
			print('Step [' + str(current_step) + '/' + total_steps + '] - Processing Search Results', end = '\r')
			start_time = time.time()
			subprocess.run([corrected_origin_path+'Python_Scripts/ProcessData/./extraction.sh', ref_result, var_result, result])
			print('Step [' + str(current_step) + '/' + total_steps + '] - Processing Search Results DONE: %.2f seconds' % round(time.time() - start_time,2))
			current_step += 1

			#Pam Analysis and creation
			print('Step [' + str(current_step) + '/' + total_steps + '] - PAM Analysis', end = '\r')
			start_time = time.time()
			print('Step [' + str(current_step) + '/' + total_steps + '] - PAM Analysis: Reference Targets', end = '\r')
			with open(result + '.semi_common_targets.minmaxdisr.txt', 'w+') as targ_result:
				subprocess.run(['awk','{real_guide=$2; gsub("-","",real_guide); print $0"\tn\tn\tn\tn\t"real_guide}', result + '.semi_common_targets.txt'], stdout = targ_result)
			
			print('\033[KStep [' + str(current_step) + '/' + total_steps + '] - PAM Analysis: Enriched Targets', end = '\r')
			with open(result + '.unique_targets.pamcreation.txt', 'w+') as targ_result:
				subprocess.run(['awk','{real_guide=$2; gsub("-","",real_guide); print $0"\tn\ty\tn\tn\t"real_guide}', result + '.unique_targets.txt'], stdout = targ_result)
			
			with open(result + '.total.txt', 'w+') as targ_result:
				subprocess.run(['cat', result + '.unique_targets.pamcreation.txt', result + '.semi_common_targets.minmaxdisr.txt'], stdout = targ_result)
			print('Step [' + str(current_step) + '/' + total_steps + '] - PAM Analysis DONE: %.2f seconds' % round(time.time() - start_time,2))
			current_step += 1
		else:
			#Pam Analysis and creation for VAR
			print('Step [' + str(current_step) + '/' + total_steps + '] - PAM Analysis', end = '\r')
			start_time = time.time()
			with open(result + '.total.txt', 'w+') as targ_result:
				subprocess.run(['awk','{real_guide=$2; gsub("-","",real_guide); print $0"\tn\tn\tn\tn\t"real_guide}', var_result], stdout = targ_result)
			print('Step [' + str(current_step) + '/' + total_steps + '] - PAM Analysis DONE: %.2f seconds' % round(time.time() - start_time,2))
			current_step += 1

		#Clustering
		print('Step [' + str(current_step) + '/' + total_steps + '] - Clustering Targets', end = '\r')
		step = 'Step [' + str(current_step) + '/' + total_steps + '] - Clustering Targets'
		start_time = time.time()
		subprocess.run([corrected_origin_path+'Python_Scripts/ProcessData/cluster.dict.py', file_to_cluster, 'no', 'True', only_cluster, guides, result + '.total', 'total', 'orderChr', step])
		#Get max mm, bulge dna and rna
		max_values = subprocess.Popen(['awk', 'BEGIN{m_m=0; b_d=0; b_r=0} NR>1 {if ($8>m_m) m_m=$8; if ($1=="DNA" && $9>b_d) b_d=$9; else if ($1=="RNA" && $9>b_r) b_r=$9} END{print m_m","b_d","b_r}', result + '.total.cluster.txt'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = max_values.communicate()
		max_values = out.decode('UTF-8').strip().split(',')	#Get number of max mm, dna and rna bulge
		max_bulge = max_values[1]
		if int(max_values[2]) > int(max_values[1]):
			max_bulge = max_values[2]
		print('Step [' + str(current_step) + '/' + total_steps + '] - Clustering Targets DONE: %.2f seconds' % round(time.time() - start_time,2))
		current_step += 1

		#Get sample and annotations and scores
		print('Step [' + str(current_step) + '/' + total_steps + '] - Samples Extraction and Annotation', end = '\r')
		step = 'Step [' + str(current_step) + '/' + total_steps + '] - Samples Extraction and Annotation'
		start_time = time.time()
		subprocess.run([corrected_origin_path+'Python_Scripts/ProcessData/' + annotator_script_file, annotationFile, result + '.total.cluster.txt', result, dict_directory, pam, max_values[0], reference_genome_dir, guides, max_values[1], max_values[2], sampleIDfile, step])
		subprocess.run(['mv', result + '.cluster.tmp.txt', result + '.total.cluster.txt'])
		print('\033[KStep [' + str(current_step) + '/' + total_steps + '] - Samples Extraction and Annotation DONE: %.2f seconds' % round(time.time() - start_time,2))
		current_step += 1

		#Summary position
		print('Step [' + str(current_step) + '/' + total_steps + '] - Creating Summaries: Summary by Position', end = '\r')
		start_time = time.time()
		subprocess.run([corrected_origin_path + 'Python_Scripts/ProcessData/summary_by_guide_position.py', result + '.total.cluster.txt', max_values[0], max_values[1], max_values[2], guides, result, 'Uniq'])

		#Summary guide
		print('\033[KStep [' + str(current_step) + '/' + total_steps + '] - Creating Summaries: Summary by Guide', end = '\r')
		subprocess.run([corrected_origin_path + 'Python_Scripts/ProcessData/summary_by_guide.py', result + '.samples.annotation.txt', max_values[0], max_values[1], max_values[2], guides, result, sum_guide_type])
		
		#Summary of samples
		print('\033[KStep [' + str(current_step) + '/' + total_steps + '] - Creating Summaries: Summary by Samples', end = '\r')
		subprocess.run([corrected_origin_path+'Python_Scripts/ProcessData/summary_by_samples.py', result + '.samples.annotation.txt', result, 'both', guides, sampleIDfile])
		print('\033[KStep [' + str(current_step) + '/' + total_steps + '] - Creating Summaries DONE: %.2f seconds' % round(time.time() - start_time,2))
		current_step += 1

		#Generate population distributions
		current_guide = 0
		print('Step [' + str(current_step) + '/' + total_steps + '] - Generating Report', end = '\r')
		start_time = time.time()
		with open(guides) as all_guides:
			for g in all_guides:
				g = g.strip()
				for i in range(int(max_bulge) + int(max_values[0]) + 1 ):
					subprocess.run([corrected_origin_path + 'Python_Scripts/ProcessData/populations_distribution.py', result + '.PopulationDistribution.txt', str(i), g])

		#Create Final file for each guide
		with open(guides) as all_guides:
			for g in all_guides:
				g = g.strip()
				current_guide += 1
				print('\033Step [' + str(current_step) + '/' + total_steps + '] - Generating Report: Guide ' + str(current_guide) + '/' + str(num_guides), end = '\r')
				#NOTE attenzione che qui viene usato shell = True
				subprocess.run(['LC_ALL=C grep ' + g + ' ' + result + '.samples.annotation.txt > ' + result + '.targets.' + g + '.txt'], shell = True)
				subprocess.run([corrected_origin_path + 'Python_Scripts/ProcessData/cluster.dict.py', result + '.targets.' + g + '.txt', 'no', 'True', 'True', guides, result + '.targets.' + g ,'total', 'addForFinal'])
				subprocess.run(['mv', result + '.targets.' + g + '.cluster.txt', result + '.targets.' + g + '.txt'])
		

		#Remove intermediate files
		file_to_remove = ['.total', '.common_targets', '.semi_common_targets','.semi_common_targets.minmaxdisr', '.unique_targets',
		'.unique_targets.pamcreation', '.addToGeneralTable','.SampleClasses', '.samples.annotation']
		for i in file_to_remove:
			removeFile(result + i + '.txt')
		print('Step [' + str(current_step) + '/' + total_steps + '] - Generating Report DONE: %.2f seconds' % round(time.time() - start_time,2))
	print('Processing Data END')
	print("Runtime: %s seconds" % (time.time() - start_time_processdata))

def graphicalInterface():
	#start graphical interface
	if ('help' in sys.argv[1:]):		#TODO completare
		print('Open the Graphical Interface for CRISPRitz.',
		"\nEXAMPLE CALL: crispritz.py graphical-interface <serverDirectory>\n",
		'\n<serverDirectory>: (Optional) directory where to start the server for the graphical interface. If not provided, the current directory will be used\n',
		'Please note that the <serverDirectory> must be structured as such:\n',
		'-annotations\n\t-<annotationsFile>: One or more text files containing the annotations in .bed format\n',
		'-assets\n',
		'-dictionaries\n\t-dictionary_name\n',
		'-genome_library\n',
		'-Genomes\n',
		'-pam\n',
		'-Results')
		sys.exit()
	try:
		launch_dir = os.path.realpath(sys.argv[2])
	except:
		launch_dir = os.path.realpath('./')
	if not os.path.isdir(launch_dir):
		print('Directory not found. Creating...')
		os.mkdir(launch_dir)
	os.chdir(launch_dir)
	check_directories = ['annotations', 'assets', 'dictionaries', 'genome_library', 'Genomes', 'pam', 'Results']
	for d in check_directories:
		if not os.path.isdir(os.path.realpath(d)):
			print(d + ' directory not found. Creating...')
			os.mkdir(os.getcwd() + '/' + d)
		
	print("LAUNCHING GRAPHICAL SERVER...")
	subprocess.run([corrected_origin_path + 'Python_Scripts/CrispritzWebApp/app_v6.py'])
	time.sleep(1)
	print("COPY THE FOLLOWING LINK INTO ANY BROWSER: http//:localhost")


#HELP FUNCTION
def callHelp():

	print("help:\n",
	# "\ncrispritz graphical-interface FUNCTION TO START LOCAL SERVER AND ACTIVATE USER GRAPHICAL INTERFACE"
	"\ncrispritz add-variants FUNCTION TO ADD VARIANTS DATA TO A FASTA GENOME",
	"\ncrispritz index-genome FUNCTION TO CREATE GENOME INDEX TO PERFORM FAST SEARCHES WITH BULGES",
	"\ncrispritz search FUNCTION TO PERFORM SEARCHES ON A GENOME (INDEXED OR PLAIN FASTA)",
	"\ncrispritz scores FUNCTION TO CALCULATE THE CFD SCORE FOR A LIST OF TARGETS",
	"\ncrispritz annotate-results FUNCTION TO ADD GENOMIC INFORMATION TO TARGETS RESULTS",
	"\ncrispritz generate-report FUNCTION TO GENERATE GRAPHICAL REPORT FOR A SPECIFIC GUIDE",
	"\ncrispritz process-data FUNCTION TO ANALYZE RESULTS TO GENERATE VARIANT ANALYSIS AND SAMPLE CLASSIFICATION",
	"\n\nADD help TO ANY FUNCTION TO VISUALIZE A BRIEF HELP PAGE (example: crispritz index-genome help)\n")

if len(sys.argv) < 2:
	callHelp()
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
elif sys.argv[1] == 'process-data':
	processData()
elif sys.argv[1] == 'graphical-interface':
	graphicalInterface()
else:
	print("ERROR! \"" + sys.argv[1] + "\" is not an allowed!")

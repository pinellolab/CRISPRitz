#!/usr/bin/env python

# Python Program to search CRISPR/Cas complex into a genome
import subprocess					# run c++ executable
import time
import os							# instructions manage directories
import shutil						# remove directory and its content
import sys							# input argv

# path where this file is located
origin_path = os.path.dirname(os.path.realpath(__file__))

corrected_origin_path = (origin_path[:-3])


def callHelp():

	print("help:",
	"\n\tcrispritz add-variants <vcfFilesDirectory> <genomeDirectory>",
	"\n\tcrispritz index-genome <name_genome> <genomeDirectory> <pamFile> -bMax <maxBulges_num>",
	"\n\tcrispritz search <genomeDirectory> <pamFile> <guidesFile> <outputFile> {-index} (flag to search with index-genome, allow searching with bulges) -mm <mm_num> [-bRNA <bRNA_num> | -bDNA <bDNA_num>] [-th <num_thread>] {-r,-p,-t} (write only off-targets results, write only profiles, write both) [-var] (to activate search with IUPAC nomenclature) [-scores <fastaGenomeDirectory> (to calculate CFD and Doench 2016 scores)]",
	"\n\tcrispritz annotate-results <guidesFile> <resultsFile> <annotationsPathFile> <outputFile>",
	"\n\tcrispritz generate-report <guide> -mm <mm_num or range mm_min-mm_max> -profile <guideProfile> -extprofile <guideExtendedProfile> -exons <exonsCountFile> -introns <intronsCountFile> -ctcf <CTCFCountFile> -dnase <DNAseCountFile> -promoters <promotersCountFile> [-gecko (to use gecko pre-computed profile)] [-sumref <summaryReferenceCountFile>][-sumenr <summaryEnrichedCountFile>]")


def indexGenome():
	if (len(sys.argv) < 7):
		print("WARNING: Too few arguments to function index-genome. Please provide:",
		"\n\n<name_genome> : Name of the genome to create",
		"\n\n<genomeDirectory> : Directory containing a genome in fasta format, need to be separated into single chromosome files.",
		"\n\n<pamFile>: Text file containing the PAM (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence",
		"\n\n-bMax <maxBulges_num> : Number of bulges allowed for the search phase")
		sys.exit()
	nameGenome = sys.argv[2]								# save name of the genome
	dirGenome = os.path.realpath(sys.argv[3])
	dirPAM = os.path.realpath(sys.argv[4])
	filePAM = open(os.path.realpath(sys.argv[4]), "r")
	listChrs = os.listdir(dirGenome)						# save list of chromosomes
	
	max_bulges = (sys.argv).index("-bMax") + 1
	max_bulges = sys.argv[max_bulges]
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

	# variant
	variant = 0
	if "-var" in sys.argv[1:]:
		variant = 1

	if os.path.isdir(dirTSTgenome):							# check if TSTgenome dir exists
		shutil.rmtree(dirTSTgenome)							# remove old TSTgenome dir
	os.makedirs(dirTSTgenome)								# build new TSTgenome dir
	os.chdir(dirTSTgenome)									# move into the TSTgenome dir

	# run buildTST
	start_time = time.time()
	for f in listChrs:
		print("Indexing:", f)
		subprocess.run([corrected_origin_path+"opt/crispritz/buildTST",
		               str(dirGenome)+"/"+str(f), str(dirPAM), str(variant), max_bulges])
	print("Finish indexing")
	print("Indexing runtime: %s seconds" % (time.time() - start_time))


def searchTST():
	
	if (len(sys.argv) < 10):
		print("WARNING: Too few arguments to function search. Please provide:",
		"\n\n<genomeDirectory> : Directory containing an index genome in .bin format, separated into single chromosome files",
		"\n\n<pamFile> : Text file containing the PAM sequence (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence" ,
		"\n\n<guidesFile> : Text file containing one or more guides (including a number of Ns equal to the length of the PAM sequence)", 
		"\n\n<outputFile> : Name of output file" ,
		"\n\n-index : Tag to activate index search",
		"\n\n-mm <mm_num> : Number of allowed mismatches", 
		"\n\n-bRNA <bRNA_num> : (Optional) Size of RNA bulges",
		"\n\n-bDNA <bDNA_num> : (Optional) Size of DNA bulges" ,
		"\n\n-th <num_thread> : (Optional) Number of threads to use",
		"\n\n{-r,-p,-t} : Output type (-r off-targets list only, -p profile only, -t everything)",
		"\n\n-var : (Optional) Tag to activate search with IUPAC nomenclature",
		"\n\n-scores <fastaGenomeDirectory> : (Optional) Tag to calculate CFD and Doench 2016 scores of the output targets. The directory containing the fasta files (.fasta or .fa) of the chromosomes is also needed ")
		sys.exit()

	nameGenome = os.path.realpath(sys.argv[2])
	PAM = os.path.realpath(sys.argv[3])								# save PAM
	fileGuide = os.path.realpath(sys.argv[4])
	nameResult = (sys.argv[5])  # name of result file
	dirTSTgenome = os.path.realpath(sys.argv[2])+"/"
	max_bulges = dirTSTgenome.split("/")[-2].split("_")[1]
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
		print("WARNING! Max available bulges (" + max_bulges + ") is smaller than input DNA bulges (" + bDNA + "). Using " + max_bulges + " bulges")
		bDNA = max_bulges
	if (int(max_bulges) < bRNA):
		print("WARNING! Max available bulges (" + max_bulges + ") is smaller than input RNA bulges (" + bRNA + "). Using " + max_bulges + " bulges")
		bRNA = max_bulges
	# read number of threads
	th = 1
	if "-th" in sys.argv[1:]:
		try:
			th = (sys.argv).index("-th") + 1
			th = int(sys.argv[th])
		except:
			print("ATTENTION! Check the mismatches option: -th <th_num> (th_num is an integer)")
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

	# run searchOnTST
	print("Search START")
	start_time = time.time()
	subprocess.run([corrected_origin_path+"opt/crispritz/searchTST", str(dirTSTgenome), str(fileGuide),
						str(mm), str(bDNA), str(bRNA), str(PAM), str(nameResult), str(r), str(th), max_bulges])
	print("Search END")
	print("Search runtime: %s seconds" % (time.time() - start_time))
	if "-scores" in sys.argv[1:]:

		PAM = os.path.realpath(sys.argv[3])
		try:
			idx_genome_fasta = (sys.argv).index("-scores") + 1
			idx_genome_fasta = os.path.realpath(sys.argv[idx_genome_fasta])
		except:
			print("Please select the directory containing the fasta files of the genome")
			sys.exit()
		
		pam_len = int(open(PAM).readline().split(" ")[1])
		if (pam_len < 0):
			pam_begin = True
		else:
			pam_begin = False
		pam_guide = len(open(PAM).readline().split(" ")[0])

		if (pam_guide != 23):
			print("WARNING: The CFD score and the Doench score can be calculated only for guides with 20nt and a 3nt PAM")
			sys.exit()
		pam_seq_check_ngg = open(PAM).readline().split(" ")[0].upper()
		if ("NGG" not in pam_seq_check_ngg):
			print("WARNING: The model used for the CFD and Doench scores are based on the NGG PAM, the scores may not be valid for other PAMs")
		target_filename = os.path.realpath(nameResult)
		subprocess.run([corrected_origin_path+'opt/crispritz/Python_Scripts/Scores/scores.py', idx_genome_fasta + "/", str(pam_guide), target_filename + ".targets.txt", str(pam_begin)])

def searchBruteForce():
	if (len(sys.argv) < 9):
		print("WARNING: Too few arguments to function search. Please provide:",
		"\n\n<genomeDirectory> : Directory containing a genome in fasta format, need to be separated into single chromosome files",
		"\n\n<pamFile> : Text file containing the PAM sequence (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence" ,
		"\n\n<guidesFile> : Text file containing one or more guides (including a number of Ns equal to the length of the PAM sequence)", 
		"\n\n<outputFile> : Name of output file" ,
		"\n\n-mm <mm_num> : Number of allowed mismatches", 
		"\n\n-th <num_thread> : (Optional) Number of threads to use",
		"\n\n{-r,-p,-t} : Output type (-r off-targets list only, -p profile only, -t everything)",
		"\n\n-var : (Optional) Tag to activate search with IUPAC nomenclature",
		"\n\n-scores <fastaGenomeDirectory> : (Optional) Tag to calculate CFD and Doench 2016 scores of the output targets. The directory containing the fasta files (.fasta or .fa) of the chromosomes is also needed ")
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
			print("ATTENTION! Check the mismatches option: -th <th_num> (th_num is an integer)")
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

	# run searchBruteForce
	print("Search START")
	start_time = time.time()
	subprocess.run([corrected_origin_path+"opt/crispritz/searchBruteForce", str(genomeDir), str(filePAM),
						str(fileGuide), str(mm), str(result), str(th), str(r), str(variant)])
	print("Search END")
	print("Search runtime: %s seconds" % (time.time() - start_time))

	if "-scores" in sys.argv[1:]:
		try:
			idx_genome_fasta = (sys.argv).index("-scores") + 1
			idx_genome_fasta = os.path.realpath(sys.argv[idx_genome_fasta])
		except:
			print("Please select the directory containing the fasta files of the genome")
			sys.exit()
		
		pam_len = int(open(filePAM).readline().split(" ")[1])
		if (pam_len < 0):
			pam_begin = True
		else:
			pam_begin = False
		 
		pam_guide = len(open(filePAM).readline().split(" ")[0])
		if (pam_guide != 23):
			print("WARNING: The CFD score and the Doench score can be calculated only for guides with 20nt and a 3nt PAM")
			sys.exit()
		pam_seq_check_ngg = open(filePAM).readline().split(" ")[0].upper()
		if ("NGG" not in pam_seq_check_ngg):
			print("WARNING: The model used for the CFD and Doench scores are based on the NGG PAM, the scores may not be valid for other PAMs")
		target_filename = os.path.realpath(result)
		subprocess.run([corrected_origin_path+'opt/crispritz/Python_Scripts/Scores/scores.py',  idx_genome_fasta + "/", str(pam_guide), target_filename + ".targets.txt", str(pam_begin)])


def annotateResults():
	if (len(sys.argv) < 6):
		print("WARNING: Too few arguments to function annotate-results. Please provide:",
		"\n\n<guidesFile> : Text file containing one or more guides, the same file used during the search phase",
		"\n\n<resultsFile> : Targets file, containing all genomic targets for the guides set",
		"\n\n<annotationsPathFile> : Text file containing the path to bed files with the annotations",
		"\n\n<outputFile> : Name of output file")
		sys.exit()
	
	guidesFile = os.path.realpath(sys.argv[2])
	resultsFile = os.path.realpath(sys.argv[3])
	annotationsFile = os.path.realpath(sys.argv[4])
	outputFile = os.path.realpath(sys.argv[5])

	inAnnotationsFile = open(annotationsFile, "r")

	exonFile = "no"
	intronFile = "no"
	ctcfFile = "no"
	promoterFile = "no"
	dnaseFile = "no"

	for line in inAnnotationsFile:
		x = line.split("\t")
		if str(x[0]) == "EXONS":
			exonFile = str(x[1]).rstrip()
		elif str(x[0]) == "INTRONS":
			intronFile = str(x[1]).rstrip()
		elif str(x[0]) == "CTCF":
			ctcfFile = str(x[1]).rstrip()
		elif str(x[0]) == "PROMOTERS":
			promoterFile = str(x[1]).rstrip()
		elif str(x[0]) == "DNASE":
			dnaseFile = str(x[1]).rstrip()
            
	print("Annotation START")
	start_time = time.time()
	subprocess.run([corrected_origin_path+'opt/crispritz/Python_Scripts/Annotator/annotator.py', str(exonFile), 
						str(intronFile), str(promoterFile), str(ctcfFile), str(dnaseFile), guidesFile, resultsFile, outputFile])
	print("Annotation END")
	print("Annotation runtime: %s seconds" % (time.time() - start_time))

def genomeEnrichment():
	if (len(sys.argv) < 4):
		print("WARNING: Too few arguments to function add-variants. Please provide:",
		"\n\n<vcfFilesDirectory> : Directory containing VCF files, need to be separated into single chromosome files (multi-sample files will be collapsed into one fake individual)", 
		"\n\n<genomeDirectory> : Directory containing a genome in fasta format, need to be separated into single chromosome files.")
		sys.exit()
	dirVCFFiles = os.path.realpath(sys.argv[2])
	dirGenome = os.path.realpath(sys.argv[3])
	listChrs = os.listdir(dirVCFFiles)

	dirParsedFiles = "./parsed_vcf_files/"

	if not (os.path.isdir(dirParsedFiles)):
		os.makedirs(dirParsedFiles)

	os.chdir(dirParsedFiles)

	print("Variants Extraction START")
	start_time = time.time()
	subprocess.run([corrected_origin_path+'opt/crispritz/Python_Scripts/Enrichment/bcf_query.sh', dirVCFFiles+'/'])
	print("Variants Extraction END")
	print("Runtime: %s seconds" % (time.time() - start_time))


	os.chdir("../")
	listChrs = os.listdir(dirParsedFiles)

	dirEnrichedGenome = "./variants_genome/"

	if (os.path.isdir(dirEnrichedGenome)):
		shutil.rmtree(dirEnrichedGenome)
		os.makedirs(dirEnrichedGenome)
	else:
		os.makedirs(dirEnrichedGenome)

	os.chdir(dirEnrichedGenome)
	os.makedirs("./SNPs_genome/")
	os.makedirs("./INDELs_genome/")

	print("Variants Adding START")
	# run buildTST
	start_time = time.time()
	for f in listChrs:
		splitf = f.split('.')
		altfile = str('../'+dirParsedFiles+splitf[0]+'.alt')
		genfile = str(dirGenome+'/'+splitf[0]+'.fa')
		print("Adding Variants to:", splitf[0])
		subprocess.run([corrected_origin_path+'opt/crispritz/Python_Scripts/Enrichment/enricher.py', altfile, genfile])

	print("Variants Adding END")
	print("Runtime: %s seconds" % (time.time() - start_time))

	os.chdir("../")
	shutil.rmtree(dirParsedFiles)

def generateReport():
	if (len(sys.argv) < 19):
		print("WARNING: Too few arguments to function generate-report. Please provide:",
		"\n\n<guide> : A guide present in the analyzed set",
		"\n\n-mm <mm_num or range mm_min-mm_max> : Number of mismatches to analyze",
		"\n\n-profile <guideProfile> : Profile file",
		"\n\n-extprofile <guideExtendedProfile> : Extended profile file",
		"\n\n-exons <exonsCountFile> : Count files for genomic annotations (ExonsCount.txt)",
		"\n\n-introns <intronsCountFile> : Count files for genomic annotations (IntronsCount.txt)",
		"\n\n-ctcf <CTCFCountFile> : Count files for genomic annotations (CTCFCount.txt)" ,
		"\n\n-dnase <DNAseCountFile> : Count files for genomic annotations (DNAseCount.txt)",
		"\n\n-promoters <promotersCountFile> : Count files for genomic annotations (PromotersCount.txt)",
		"\n\n-gecko : (Optional)Tag to activate gecko dataset comparison",
		"\n\n-sumref <summaryReferenceCountFile> : (Optional) Create a barplot with Reference genome",
		"\n\n-sumenr <summaryEnrichedCountFile> : (Optional) Create a barplot with Enriched genome")
		sys.exit()

	guidesFile = str(sys.argv[2])

	mm = 0
	if "-mm" in sys.argv[1:]:
		try:
			mm = (sys.argv).index("-mm") + 1
			mm = (sys.argv[mm])
		except:
			print("ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is a number)")
			sys.exit()

	profileFile = "no"
	if "-profile" in sys.argv[1:]:
		profileFile = (sys.argv).index("-profile") + 1
		profileFile = os.path.realpath(sys.argv[profileFile])

	extProfileFile = "no"
	if "-extprofile" in sys.argv[1:]:
		extProfileFile = (sys.argv).index("-extprofile") + 1
		extProfileFile = os.path.realpath(sys.argv[extProfileFile])

	exonFile = "no"
	if "-exons" in sys.argv[1:]:
		exonFile = (sys.argv).index("-exons") + 1
		exonFile = os.path.realpath(sys.argv[exonFile])

	intronFile = "no"
	if "-introns" in sys.argv[1:]:
		intronFile = (sys.argv).index("-introns") + 1
		intronFile = os.path.realpath(sys.argv[intronFile])

	promoterFile = "no"
	if "-promoters" in sys.argv[1:]:
		promoterFile = (sys.argv).index("-promoters") + 1
		promoterFile = os.path.realpath(sys.argv[promoterFile])

	ctcfFile = "no"
	if "-ctcf" in sys.argv[1:]:
		ctcfFile = (sys.argv).index("-ctcf") + 1
		ctcfFile = os.path.realpath(sys.argv[ctcfFile])

	dnaseFile = "no"
	if "-dnase" in sys.argv[1:]:
		dnaseFile = (sys.argv).index("-dnase") + 1
		dnaseFile = os.path.realpath(sys.argv[dnaseFile])

	summaryOne = "no"
	if "-sumref" in sys.argv[1:]:
		summaryOne = (sys.argv).index("-sumref") + 1
		summaryOne = os.path.realpath(sys.argv[summaryOne])

	summaryTwo = "no"
	if "-sumenr" in sys.argv[1:]:
		summaryTwo = (sys.argv).index("-sumenr") + 1
		summaryTwo = os.path.realpath(sys.argv[summaryTwo])

	geckoProfile = "no"
	geckoExonsCount = "no"
	geckoIntronsCount = "no"
	geckoPromotersCount = "no"
	geckoDNAseCount = "no"
	geckoCTCFCount = "no"

	if "-gecko" in sys.argv[1:]:
		geckoProfile = corrected_origin_path + \
				'opt/crispritz/Python_Scripts/Plot/gecko/gecko.reference.profile.xls'
		geckoExonsCount = corrected_origin_path + \
				'opt/crispritz/Python_Scripts/Plot/gecko/gecko.Exons.Count.txt'
		geckoIntronsCount = corrected_origin_path + \
				'opt/crispritz/Python_Scripts/Plot/gecko/gecko.Introns.Count.txt'
		geckoPromotersCount = corrected_origin_path + \
				'opt/crispritz/Python_Scripts/Plot/gecko/gecko.Promoters.Count.txt'
		geckoDNAseCount = corrected_origin_path + \
				'opt/crispritz/Python_Scripts/Plot/gecko/gecko.DNAse.Count.txt'
		geckoCTCFCount = corrected_origin_path + \
				'opt/crispritz/Python_Scripts/Plot/gecko/gecko.CTCF.Count.txt'

	subprocess.run([corrected_origin_path + 'opt/crispritz/Python_Scripts/Plot/radar_chart.py', str(profileFile), str(extProfileFile), str(exonFile), 
						str(intronFile), str(promoterFile), str(dnaseFile),str(ctcfFile), guidesFile, str(mm), str(geckoProfile), 
						str(geckoExonsCount), str(geckoIntronsCount), str(geckoPromotersCount), str(geckoDNAseCount), str(geckoCTCFCount),
						str(summaryOne), str(summaryTwo)])


if len(sys.argv) < 2:
	callHelp()
elif sys.argv[1] == "index-genome":
	indexGenome()
elif sys.argv[1] == "search" and "-index" in sys.argv[1:]:
	searchTST()
elif sys.argv[1] == "search":
	searchBruteForce()
elif sys.argv[1] == "add-variants":
	genomeEnrichment()
elif sys.argv[1] == "annotate-results":
	annotateResults()
elif sys.argv[1] == "generate-report":
	generateReport()
else:
	print("ERROR! \"" + sys.argv[1] + "\" is not an allowed!")

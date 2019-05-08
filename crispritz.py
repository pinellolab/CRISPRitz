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
	"\n\tcrispritz index-genome <name_genome> <genomeDirectory> <pamFile> ",
	"\n\tcrispritz search <genomeDirectory> <pamFile> <guidesFile> <outputFile> {-index} (flag to search with index-genome, allow searching with bulges) -mm <mm_num> [-bRNA <bRNA_num> | -bDNA <bDNA_num>] [-th <num_thread>] {-r,-p,-t} (write only off-targets results, write only profiles, write both) [-var] (to activate search with IUPAC nomenclature)",
	"\n\tcrispritz annotate-results <guidesFile> <resultsFile> <outputFile> -exons <exonsbedFile> -introns <intronsbedFile> -ctcf <ctcfbedFile> -dnase <dnasebedFile> -promoters <promotersbedFile>",
	"\n\tcrispritz generate-report <guide> -mm <mm_num or range mm_min-mm_max> -profile <guideProfile> -extprofile <guideExtendedProfile> -exons <exonsCountFile> -introns <intronsCountFile> -ctcf <CTCFCountFile> -dnase <DNAseCountFile> -promoters <promotersCountFile> [-gecko (to use gecko pre-computed profile)] [-sumref <summaryReferenceCountFile>][-sumenr <summaryEnrichedCountFile>]")

def indexGenome():

	nameGenome = sys.argv[2]								# save name of the genome
	dirGenome = os.path.realpath(sys.argv[3])
	dirPAM = os.path.realpath(sys.argv[4])
	filePAM = open(os.path.realpath(sys.argv[4]), "r")
	listChrs = os.listdir(dirGenome)						# save list of chromosomes

	# retrive PAM
	PAM = filePAM.read()
	PAM_size = int(PAM.split()[1])
	if(PAM_size < 0):
		PAM_size = PAM_size * -1
		PAM = PAM.split()[0][0:PAM_size]
	else:
		PAM = PAM.split()[0][-PAM_size:]

	TSTgenome = PAM + "_" + nameGenome					# name of the genome in TST format
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
		subprocess.run([corrected_origin_path+"opt/crispritz/buildTST", str(dirGenome)+"/"+str(f), str(dirPAM), str(variant)])
	print("Finish indexing")
	print("Indexing runtime: %s seconds" % (time.time() - start_time))

def searchTST():

	nameGenome = sys.argv[2]
	PAM = sys.argv[3]								# save PAM
	fileGuide = os.path.realpath(sys.argv[4])
	nameResult = (sys.argv[5])  # name of result file
	dirTSTgenome = os.path.realpath(sys.argv[2])+"/"

	if not os.path.isdir(dirTSTgenome):				# check if TSTgenome dir exists
		print("ATTENTION! You have to generate the index of \"" + nameGenome + "\" with \"" + PAM + "\", before the search using index!")
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
						str(mm), str(bDNA), str(bRNA), str(PAM), str(nameResult), str(r), str(th)])
	print("Search END")
	print("Search runtime: %s seconds" % (time.time() - start_time))

def searchBruteForce():

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

def annotateResults():

	guidesFile = os.path.realpath(sys.argv[2])
	resultsFile = os.path.realpath(sys.argv[3])
	outputFile = os.path.realpath(sys.argv[4])

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

	print("Annotation START")
	start_time = time.time()
	subprocess.run([corrected_origin_path+'opt/crispritz/Python_Scripts/Annotator/annotator.py', str(exonFile), 
						str(intronFile), str(promoterFile), str(ctcfFile), str(dnaseFile), guidesFile, resultsFile, outputFile])
	print("Annotation END")
	print("Annotation runtime: %s seconds" % (time.time() - start_time))

def genomeEnrichment():
	
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

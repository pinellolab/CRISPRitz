#!/usr/bin/python3
# Python Program to search CRISPR/Cas complex into a genome
import subprocess					# run c++ executable
import time
import os							# instructions manage directories
import shutil						# remove directory and its content
import sys							# input argv


# options:
mismatches = ["--mismatches", "-mm"]				# mismatches
bulgeRNA = ["--bulgeRNA", "-bRNA"]					# bulge RNA
bulgeDNA = ["--bulgeDNA", "-bDNA"]					# bulge DNA


if len(sys.argv) < 2:
	print("help:",
	"\n\tcrispritz index-gen <name_genome> <genomeDirectory> <pamFile> ",
	"\n\tcrispritz index-search <name_genome> <PAM> <guideFile> (-mm <mm_num> | -bRNA <bRNA_num> | -bDNA <bDNA_num>)",
	"\n\tcrispritz search <genomeDirectory> <pamFile> <guideFile> <resultFile> -mm <mm_num> [-th <num_thread>] {-r,-p,-t} (write only off-targets results,write only profiles, write both)",
	"\n\tcrispritz enrichment <vcf_parsedDirectory> <genomeDirectory>",
	"\n\tcrispritz vcf-parsing <vcfFilesDirectory>",
	"\n\tcrispritz annotation <guidesFile> <resultsFile> <outputFile> -exons <exonsbedFile> -introns <intronsbedFile> -ctcf <ctcfbedFile> -dnase <dnasebedFile> -promoters <promotersbedFile>",
	"\n\tcrispritz guide-report <guide> -mm <mm_num> [-profile <guideProfile> -extprofile <guideExtendedProfile> -exons <exonsCountFile> -introns <intronsCountFile> -ctcf <CTCFCountFile> -dnase <DNAseCountFile> -promoters <promotersCountFile>][-gecko (to use gecko pre-computed profile), -random1milion (to use random 1 milion pre-computed profile)]")
	
	
else:
	# -------------------- Index generation --------------------
	if sys.argv[1] == "index-gen":
		nameGenome = sys.argv[2]								# save name of the genome
		dirGenome = os.path.realpath(sys.argv[3])				# save dir of the genome (with absolute path)
		filePAM = open(os.path.realpath(sys.argv[4]), "r")		# open file PAM (with absolute path)
		listChrs = os.listdir(dirGenome)						# save list of chromosomes
		
		# retrive PAM
		PAM = filePAM.read()
		PAM_size = int(PAM.split()[1])
		PAM = PAM.split()[0][-PAM_size:]
		
		TSTgenome = PAM + "_" + nameGenome					# name of the genome in TST format
		dirTSTgenome = "./genome_library/" + TSTgenome		# dir of the genome in TST format
		
		print (TSTgenome, "Indexing generation:")
		
		# compile buildTST
		if (os.path.isfile("buildTST")) == False:
			print ("Compiling...")
			subprocess.run(["g++", "-O3", "sourceCode/CRISPR-Cas-Tree/main.cpp", "-o", "buildTST", "-fopenmp", "-std=c++11"])
		
		if os.path.isdir(dirTSTgenome):							# check if TSTgenome dir exists
			shutil.rmtree(dirTSTgenome)							# remove old TSTgenome dir
		os.makedirs(dirTSTgenome)								# build new TSTgenome dir
		os.chdir(dirTSTgenome)									# move into the TSTgenome dir
		
		# run buildTST
		start_time = time.time()
		for f in listChrs:
			print ("Indexing:", f)
			subprocess.run(["./../../buildTST", str(dirGenome)+"/"+str(f), PAM])
		print ("Finish indexing")
		print("Indexing runtime: %s seconds" % (time.time() - start_time))
		
	
	
	# -------------------- Indexed search --------------------
	elif sys.argv[1] == "index-search":
		
		nameGenome = sys.argv[2]						# save name of the genome
		PAM = sys.argv[3]								# save PAM
		fileGuide = os.path.realpath(sys.argv[4])		# open file Guide (with absolute path)
		TSTgenome = PAM + "_" + nameGenome				# name of the genome in TST format
		dirTSTgenome = "./genome_library/" + TSTgenome	# dir of the genome in TST format
		
		if not os.path.isdir(dirTSTgenome):				# check if TSTgenome dir exists
			print ("ATTENTION! You have to generate the index of \"" + nameGenome + "\" with \"" + PAM + "\", before the search using index!")
			sys.exit()
		
		listChrs = os.listdir(dirTSTgenome)				# save list of chromosomes
		
		print (TSTgenome, "Search using indexing:")
		
		# compile searchOnTST
		if (os.path.isfile("searchOnTST")) == False:
			print ("Compiling...")
			subprocess.run(["g++", "-O3", "sourceCode/CRISPR-Cas-Tree/searchOnTST.cpp", "-o", "searchOnTST", "-fopenmp", "-std=c++11"]) 
		
		# read number of mismatches
		mm = 0
		if "-mm" in sys.argv[1:]:
			try:
				mm = (sys.argv).index("-mm") + 1
				mm = int(sys.argv[mm])
			except:
				print ("ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is a number)")
				sys.exit()
		
		# read number of bulge RNA
		bRNA = 0
		if "-bRNA" in sys.argv[1:]:
			try:
				bRNA = (sys.argv).index("-bRNA") + 1
				bRNA = int(sys.argv[bRNA])
			except:
				print ("ATTENTION! Check the bulge RNA option: -bRNA <bRNA_num> (bRNA_num is a number)")
				sys.exit()
				
		# read number of bulge DNA
		bDNA = 0
		if "-bDNA" in sys.argv[1:]:
			try:
				bDNA = (sys.argv).index("-bDNA") + 1
				bDNA = int(sys.argv[bDNA])
			except:
				print ("ATTENTION! Check the bulge DNA option: -bDNA <bDNA_num> (bDNA_num is a number)")
				sys.exit()
		
		# write the header of results file
		text_file = open("result.txt", "w")
		text_file.write("#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge Size\n")
		text_file.close()
		
		# run searchOnTST
		start_time = time.time()
		for f in listChrs:
			print ("Search on", f)
			subprocess.run(["./searchOnTST", str(dirTSTgenome)+"/"+str(f), str(fileGuide), str(mm), str(bDNA), str(bRNA)])
		print ("Finish Searching")
		print("Searching runtime: %s seconds" % (time.time() - start_time))
				
	

	
	# -------------------- Normal search --------------------
	elif sys.argv[1] == "search":
		
		genomeDir = os.path.realpath(sys.argv[2])+"/"	# genome directory (with absolute path)
		filePAM = os.path.realpath(sys.argv[3])		# file PAM (with absolute path)
		fileGuide = os.path.realpath(sys.argv[4])	# file Guide (with absolute path)
		#numThread = sys.argv[6]						# number of allowed threads
		result = sys.argv[5]						#name of results
		
		# compile crisprofiler
		if (os.path.isfile("crisProfiler")) == False:
			print ("Compiling...")
			source = "sourceCode/CRISPRofiler/"
			subprocess.run(["g++", "-std=c++11", "-O3", "-fopenmp", str(source)+"main.cpp", str(source)+"analysis.cpp", str(source)+"profiling.cpp", str(source)+"guide_searching.cpp", str(source)+"pam_searching.cpp", str(source)+"pre_computation.cpp", str(source)+"reading.cpp", "-o", "crisProfiler"])
		
		# read number of mismatches
		mm = 0
		if "-mm" in sys.argv[1:]:
			try:
				mm = (sys.argv).index("-mm") + 1
				mm = int(sys.argv[mm])
			except:
				print ("ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is an integer)")
				sys.exit()

		# read number of mismatches
		th = 1000
		if "-th" in sys.argv[1:]:
			try:
				th = (sys.argv).index("-th") + 1
				th = int(sys.argv[th])
			except:
				print ("ATTENTION! Check the mismatches option: -th <th_num> (th_num is an integer)")
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
		print ("Start Searching")
		start_time = time.time()
		subprocess.run(["./crisProfiler", str(genomeDir), str(filePAM), str(fileGuide), str(mm), str(result), str(th), str(r)])
		print ("Finish Searching")
		print("Searching runtime: %s seconds" % (time.time() - start_time))
	
	# -------------------- Genome Enrichment --------------------
	elif sys.argv[1] == "enrichment":

		dirEnrichedFiles = os.path.realpath(sys.argv[2])								# save name of the genome
		dirGenome = os.path.realpath(sys.argv[3])				# save dir of the genome (with absolute path)
		listChrs = os.listdir(dirEnrichedFiles)						# save list of chromosomes
		
		# retrive PAM	
		dirEnrichedGenome = "./enriched_genome/"		# dir of the genome in TST format
		
		if (os.path.isdir(dirEnrichedGenome)) == False:	
			os.makedirs(dirEnrichedGenome)						# check if TSTgenome dir exists
																#	shutil.rmtree(dirEnrichedGenome)
																# remove old TSTgenome di
		os.chdir(dirEnrichedGenome)

		print ("Start Enrichment")
		# run buildTST
		start_time = time.time()
		for f in listChrs:
			splitf = f.split('.')
			altfile = str(dirEnrichedFiles+'/'+splitf[0]+'.alt')
			genfile = str(dirGenome+'/'+splitf[0]+'.fa')
			print ("Enriching:", splitf[0])
			#print ('Alt file', altfile)
			#subprocess.run(['../sourceCode/Python_Scripts/Enrichment/./enricher.sh', altfile, genfile])
			subprocess.run(['../sourceCode/Python_Scripts/Enrichment/./enricher.py', altfile, genfile])


		print ("Finish Enrichment")
		print("Enrichment runtime: %s seconds" % (time.time() - start_time))
	
	
	# -------------------- VCF File Parsing --------------------
	elif sys.argv[1] == "vcf-parsing":

		dirVCFFiles = os.path.realpath(sys.argv[2])								# save name of the genome
		#dirGenome = os.path.realpath(sys.argv[3])				# save dir of the genome (with absolute path)
		listChrs = os.listdir(dirVCFFiles)						# save list of chromosomes
		
		# retrive PAM	
		dirEnrichedGenome = "./Parsed_VCF_Files/"		# dir of the genome in TST format
		
		if (os.path.isdir(dirEnrichedGenome)) == False:	
			os.makedirs(dirEnrichedGenome)						# check if TSTgenome dir exists
																#	shutil.rmtree(dirEnrichedGenome)
																# remove old TSTgenome di
		os.chdir(dirEnrichedGenome)

		print ("Start Parsing")
		# run buildTST
		start_time = time.time()
		#for f in listChrs:
			#splitf = f.split('.')
			#altfile = str(dirEnrichedFiles+'/'+splitf[0]+'.alt')
			#genfile = str(dirGenome+'/'+splitf[0]+'.fa')
			#print ("Parsing:", f)
			#print ('Alt file', altfile)
		subprocess.run(['../sourceCode/Python_Scripts/Enrichment/./bcf_query.sh', dirVCFFiles+'/'])

		print ("Finish Parsing")
		print("Parsing runtime: %s seconds" % (time.time() - start_time))

	
	# -------------------- Annotator --------------------
	elif sys.argv[1] == "annotation":

		#exonFile = os.path.realpath(sys.argv[2])	# genome directory (with absolute path)
		#intronFile = os.path.realpath(sys.argv[3])		# file PAM (with absolute path)
		guidesFile = os.path.realpath(sys.argv[2])
		resultsFile = os.path.realpath(sys.argv[3])	
		outputFile = os.path.realpath(sys.argv[4])		# save dir of the genome (with absolute path)
		
		# retrive PAM	
		#dirEnrichedGenome = "./Parsed_VCF_Files/"		# dir of the genome in TST format
		
	#	if (os.path.isdir(dirEnrichedGenome)) == False:	
	#		os.makedirs(dirEnrichedGenome)						# check if TSTgenome dir exists
																#	shutil.rmtree(dirEnrichedGenome)
																# remove old TSTgenome di
	#	os.chdir(dirEnrichedGenome)

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


		print ("Start Annotating")
		# run buildTST
		start_time = time.time()
		#for f in listChrs:
			#splitf = f.split('.')
			#altfile = str(dirEnrichedFiles+'/'+splitf[0]+'.alt')
			#genfile = str(dirGenome+'/'+splitf[0]+'.fa')
			#print ("Parsing:", f)
			#print ('Alt file', altfile)
		#subprocess.run(['sourceCode/Python_Scripts/Annotator/./annotator.sh', exonFile, intronFile, resultsFile, outputFile])
		subprocess.run(['sourceCode/Python_Scripts/Annotator/./annotator.py', str(exonFile), str(intronFile),str(promoterFile),str(ctcfFile),str(dnaseFile), guidesFile, resultsFile, outputFile])

		print ("Finish Annotating")
		print("Annotation runtime: %s seconds" % (time.time() - start_time))


	# -------------------- Guide Report --------------------
	elif sys.argv[1] == "guide-report":

		guidesFile = str(sys.argv[2])

		mm = 0
		if "-mm" in sys.argv[1:]:
			try:
				mm = (sys.argv).index("-mm") + 1
				mm = (sys.argv[mm])
			except:
				print ("ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is a number)")
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

		geckoProfile = "no"
		geckoExonsCount = "no"
		geckoIntronsCount = "no"
		geckoPromotersCount = "no"
		geckoDNAseCount = "no"
		geckoCTCFCount = "no"
		if "-gecko" in sys.argv[1:]:
			geckoProfile = 'sourceCode/Python_Scripts/Plot/gecko/gecko.reference.profile.xls'
			geckoExonsCount = 'sourceCode/Python_Scripts/Plot/gecko/gecko.Exons.Count.txt'
			geckoIntronsCount = 'sourceCode/Python_Scripts/Plot/gecko/gecko.Introns.Count.txt'
			geckoPromotersCount = 'sourceCode/Python_Scripts/Plot/gecko/gecko.Promoters.Count.txt'
			geckoDNAseCount = 'sourceCode/Python_Scripts/Plot/gecko/gecko.DNAse.Count.txt'
			geckoCTCFCount = 'sourceCode/Python_Scripts/Plot/gecko/gecko.CTCF.Count.txt'

		subprocess.run(['sourceCode/Python_Scripts/Plot/./radar_chart.py', str(profileFile), str(extProfileFile),str(exonFile),str(intronFile),str(promoterFile),str(dnaseFile),str(ctcfFile),guidesFile,str(mm),str(geckoProfile),str(geckoExonsCount),str(geckoIntronsCount),str(geckoPromotersCount),str(geckoDNAseCount),str(geckoCTCFCount)])

	
	# -------------------- Error --------------------
	else:
		print ("ERROR! \"" + sys.argv[1] + "\" is not an allowed!")

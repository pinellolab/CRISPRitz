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


if len(sys.argv) == 1:
	print("help: \n\tcrispritz index-gen <name_genome> <path/genomeDirectory> <path/pamFile> ",
	"\n\tcrispritz index-search <name_genome> <PAM> <path/guideFile> (-mm <mm_num> | -bRNA <bRNA_num> | -bDNA <bDNA_num>)",
	"\n\tcrispritz search <path/genomeDirectory> <path/pamFile> <path/guideFile> <num_thread> (-mm <mm_num> | -r) (-r: write the results)")
	
	
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
			subprocess.run(["./searchOnTST", str(dirTSTgenome)+"/"+str(f), str(fileGuide), str(mm), str(bRNA), str(bDNA)])
		print ("Finish Searching")
		print("Searching runtime: %s seconds" % (time.time() - start_time))
				
	

	
	# -------------------- Normal search --------------------
	elif sys.argv[1] == "search":
		
		genomeDir = os.path.realpath(sys.argv[2])+"/"	# genome directory (with absolute path)
		filePAM = os.path.realpath(sys.argv[3])		# file PAM (with absolute path)
		fileGuide = os.path.realpath(sys.argv[4])	# file Guide (with absolute path)
		numThread = sys.argv[5]						# number of allowed threads
		output = os.getcwd() + "/result"		# output 
		
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
				print ("ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is a number)")
				sys.exit()
		
		# don't write the resutls
		r = ""
		if "-r" in sys.argv[1:]:
			r = "r"
		
		# run searchOnTST
		print ("Start Searching")
		start_time = time.time()
		subprocess.run(["./crisProfiler", str(genomeDir), str(filePAM), str(fileGuide), str(mm), str(output), str(numThread), str(r)])
		print ("Finish Searching")
		print("Searching runtime: %s seconds" % (time.time() - start_time))
	
	
	
	
	
	# -------------------- Error --------------------
	else:
		print ("ERROR! \"" + sys.argv[1] + "\" is not an allowed!")
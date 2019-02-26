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
      "\n\tcrispritz add-variants <vcfFilesDirectory> <genomeDirectory>",
      "\n\tcrispritz index-genome <name_genome> <genomeDirectory> <pamFile> ",
      "\n\tcrispritz search <genomeDirectory> <pamFile> <guidesFile> <outputFile> {-index} (flag to search with index-genome, allow searching with bulges) -mm <mm_num> [-bRNA <bRNA_num> | -bDNA <bDNA_num>] [-th <num_thread>] {-r,-p,-t} (write only off-targets results, write only profiles, write both)",
      "\n\tcrispritz annotate-results <guidesFile> <resultsFile> <outputFile> -exons <exonsbedFile> -introns <intronsbedFile> -ctcf <ctcfbedFile> -dnase <dnasebedFile> -promoters <promotersbedFile>",
      "\n\tcrispritz generate-report <guide> -mm <mm_num or range mm_min-mm_max> -profile <guideProfile> -extprofile <guideExtendedProfile> -exons <exonsCountFile> -introns <intronsCountFile> -ctcf <CTCFCountFile> -dnase <DNAseCountFile> -promoters <promotersCountFile> [-gecko (to use gecko pre-computed profile)] [-sumone <summaryReferenceCountFile>][-sumtwo <summaryEnrichedCountFile>]")


else:
    # -------------------- Index generation --------------------
    if sys.argv[1] == "index-genome":
        nameGenome = sys.argv[2]								# save name of the genome
        # save dir of the genome (with absolute path)
        dirGenome = os.path.realpath(sys.argv[3])
        # open file PAM (with absolute path)
        dirPAM = os.path.realpath(sys.argv[4])
        filePAM = open(os.path.realpath(sys.argv[4]), "r")
        listChrs = os.listdir(dirGenome)						# save list of chromosomes

        # retrive PAM
        PAM = filePAM.read()
        PAM_size = int(PAM.split()[1])
        PAM = PAM.split()[0][-PAM_size:]

        TSTgenome = PAM + "_" + nameGenome					# name of the genome in TST format
        dirTSTgenome = "./genome_library/" + TSTgenome		# dir of the genome in TST format

        print(TSTgenome, "Indexing generation:")

        # compile buildTST
        if (os.path.isfile("buildTST")) == False:
            print("Compiling...")
            subprocess.run(["g++", "-O3", "sourceCode/CRISPR-Cas-Tree/main.cpp","-o", "buildTST", "-fopenmp", "-std=c++11"])

        if os.path.isdir(dirTSTgenome):							# check if TSTgenome dir exists
            shutil.rmtree(dirTSTgenome)							# remove old TSTgenome dir
        os.makedirs(dirTSTgenome)								# build new TSTgenome dir
        os.chdir(dirTSTgenome)									# move into the TSTgenome dir

        # run buildTST
        start_time = time.time()
        for f in listChrs:
            print("Indexing:", f)
            subprocess.run(["./../../buildTST", str(dirGenome)+"/"+str(f), str(dirPAM)])
        print("Finish indexing")
        print("Indexing runtime: %s seconds" % (time.time() - start_time))

    # -------------------- Indexed search --------------------
    elif sys.argv[1] == "search" and "-index" in sys.argv[1:]:

        #nameGenome = sys.argv[2]						# save name of the genome
        PAM = sys.argv[3]								# save PAM
        # open file Guide (with absolute path)
        fileGuide = os.path.realpath(sys.argv[4])
        #print(fileGuide)
        #TSTgenome = PAM + "_" + nameGenome				# name of the genome in TST format
        #dirTSTgenome = "./genome_library/" + TSTgenome  # dir of the genome in TST format
        nameResult = (sys.argv[5])      #name of result file
        dirTSTgenome = os.path.realpath(sys.argv[2])+"/"

        if not os.path.isdir(dirTSTgenome):				# check if TSTgenome dir exists
            print("ATTENTION! You have to generate the index of \"" + nameGenome +
                  "\" with \"" + PAM + "\", before the search using index!")
            sys.exit()

        #listChrs = os.listdir(dirTSTgenome)				# save list of chromosomes

        #print(TSTgenome, "Search using indexing:")

        # compile searchOnTST
        if (os.path.isfile("searchOnTST")) == False:
            print("Compiling...")
            subprocess.run(["g++", "-O3", "sourceCode/CRISPR-Cas-Tree/searchTest.cpp", "sourceCode/CRISPR-Cas-Tree/detailedOutput.cpp", "sourceCode/CRISPR-Cas-Tree/convert.cpp", "-I", "sourceCode/CRISPR-Cas-Tree/include", "-o", "searchOnTST", "-fopenmp", "-std=c++11", "-g"]) 


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

        #read number of threads
        th = 1
        if "-th" in sys.argv[1:]:
            try:
                th = (sys.argv).index("-th") + 1
                th = int(sys.argv[th])
            except:
                print(
                    "ATTENTION! Check the mismatches option: -th <th_num> (th_num is an integer)")
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

        # write the header of results file
        #text_file = open("result.txt", "w")
        #text_file.write("#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge Size\n")
        #text_file.close()

        # run searchOnTST
        start_time = time.time()

        print("Search on index-genome")
        subprocess.run(["./searchOnTST", str(dirTSTgenome), str(fileGuide), str(mm), str(bDNA), str(bRNA), str(PAM), str(nameResult), str(r), str(th)]) 
        print("Finish Searching")
        print("Searching runtime: %s seconds" % (time.time() - start_time))

    # -------------------- Normal search --------------------
    elif sys.argv[1] == "search":

        # genome directory (with absolute path)
        genomeDir = os.path.realpath(sys.argv[2])+"/"
        # file PAM (with absolute path)
        filePAM = os.path.realpath(sys.argv[3])
        # file Guide (with absolute path)
        fileGuide = os.path.realpath(sys.argv[4])
        # numThread = sys.argv[6]						# number of allowed threads
        result = sys.argv[5]  # name of results

        # compile crisprofiler
        if (os.path.isfile("searchBruteForce")) == False:
            print("Compiling...")
            source = "sourceCode/CRISPRofiler/"
            subprocess.run(["g++", "-std=c++11", "-O3", "-fopenmp", str(source)+"main.cpp", str(source)+"analysis.cpp", str(source)+"profiling.cpp", str(
                source)+"guide_searching.cpp", str(source)+"pam_searching.cpp", str(source)+"pre_computation.cpp", str(source)+"reading.cpp", "-o", "searchBruteForce"])

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
        th = 10000000
        if "-th" in sys.argv[1:]:
            try:
                th = (sys.argv).index("-th") + 1
                th = int(sys.argv[th])
            except:
                print(
                    "ATTENTION! Check the mismatches option: -th <th_num> (th_num is an integer)")
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
        print("Start Searching")
        start_time = time.time()
        subprocess.run(["./searchBruteForce", str(genomeDir), str(filePAM),str(fileGuide), str(mm), str(result), str(th), str(r)])
        print("Finish Searching")
        print("Searching runtime: %s seconds" % (time.time() - start_time))

        # -------------------- VCF File Parsing --------------------
    elif sys.argv[1] == "add-variants":

        # save name of the genome
        dirVCFFiles = os.path.realpath(sys.argv[2])
        # save dir of the genome (with absolute path)
        dirGenome = os.path.realpath(sys.argv[3])
        listChrs = os.listdir(dirVCFFiles)						# save list of chromosomes

        # retrive PAM
        dirParsedFiles = "./parsed_vcf_files/"		# dir of the genome in TST format

        if (os.path.isdir(dirParsedFiles)) == False:
            os.makedirs(dirParsedFiles)						# check if TSTgenome dir exists
            #	shutil.rmtree(dirEnrichedGenome)
            # remove old TSTgenome di
        os.chdir(dirParsedFiles)

        print("Start Extracting Variants")
        # run buildTST
        start_time = time.time()
        # for f in listChrs:
        # splitf = f.split('.')
        # altfile = str(dirEnrichedFiles+'/'+splitf[0]+'.alt')
        # genfile = str(dirGenome+'/'+splitf[0]+'.fa')
        # print ("Parsing:", f)
        # print ('Alt file', altfile)
        subprocess.run(['../sourceCode/Python_Scripts/Enrichment/./bcf_query.sh', dirVCFFiles+'/'])

        print("Finish Extraction of Variants")
        print("Runtime: %s seconds" % (time.time() - start_time))
    # -------------------- Genome Enrichment --------------------

        # dirEnrichedFiles = os.path.realpath(sys.argv[2])								# save name of the genome
        # save dir of the genome (with absolute path)
        # dirGenome = os.path.realpath(sys.argv[3])
        os.chdir("../")
        listChrs = os.listdir(dirParsedFiles)						# save list of chromosomes

        # retrive PAM
        dirEnrichedGenome = "./variants_genome/"		# dir of the genome in TST format

        if (os.path.isdir(dirEnrichedGenome)) == False:
            os.makedirs(dirEnrichedGenome)						# check if TSTgenome dir exists
        os.chdir(dirEnrichedGenome)
        os.makedirs("./SNPs_genome/")
        os.makedirs("./INDELs_genome/")

        print("Start Adding Variants")
        # run buildTST
        start_time = time.time()
        for f in listChrs:
            splitf = f.split('.')
            altfile = str('../'+dirParsedFiles+splitf[0]+'.alt')
            genfile = str(dirGenome+'/'+splitf[0]+'.fa')
            print("Adding Variants to:", splitf[0])
            # print ('Alt file', altfile)
            # subprocess.run(['../sourceCode/Python_Scripts/Enrichment/./enricher.sh', altfile, genfile])
            subprocess.run(['../sourceCode/Python_Scripts/Enrichment/enricher.py', altfile, genfile])

        print("Finish Adding of Variants")
        print("Runtime: %s seconds" % (time.time() - start_time))

        os.chdir("../")
        shutil.rmtree(dirParsedFiles)

    # -------------------- Annotator --------------------
    elif sys.argv[1] == "annotate-results":

        # exonFile = os.path.realpath(sys.argv[2])	# genome directory (with absolute path)
        # intronFile = os.path.realpath(sys.argv[3])		# file PAM (with absolute path)
        guidesFile = os.path.realpath(sys.argv[2])
        resultsFile = os.path.realpath(sys.argv[3])
        # save dir of the genome (with absolute path)
        outputFile = os.path.realpath(sys.argv[4])

        # retrive PAM
        # dirEnrichedGenome = "./Parsed_VCF_Files/"		# dir of the genome in TST format

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

        print("Start Annotating")
        # run buildTST
        start_time = time.time()
        # for f in listChrs:
        # splitf = f.split('.')
        # altfile = str(dirEnrichedFiles+'/'+splitf[0]+'.alt')
        # genfile = str(dirGenome+'/'+splitf[0]+'.fa')
        # print ("Parsing:", f)
        # print ('Alt file', altfile)
        # subprocess.run(['sourceCode/Python_Scripts/Annotator/./annotator.sh', exonFile, intronFile, resultsFile, outputFile])
        subprocess.run(['sourceCode/Python_Scripts/Annotator/./annotator.py', str(exonFile), str(intronFile),
                        str(promoterFile), str(ctcfFile), str(dnaseFile), guidesFile, resultsFile, outputFile])

        print("Finish Annotating")
        print("Annotation runtime: %s seconds" % (time.time() - start_time))

    # -------------------- Guide Report --------------------
    elif sys.argv[1] == "generate-report":
      guidesFile = str(sys.argv[2])

      mm = 0
      if "-mm" in sys.argv[1:]:
         try:
               mm = (sys.argv).index("-mm") + 1
               mm = (sys.argv[mm])
         except:
               print(
                  "ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is a number)")
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
      if "-sumone" in sys.argv[1:]:
         summaryOne = (sys.argv).index("-sumone") + 1
         summaryOne = os.path.realpath(sys.argv[summaryOne])
      
      summaryTwo = "no"
      if "-sumtwo" in sys.argv[1:]:
         summaryTwo = (sys.argv).index("-sumtwo") + 1
         summaryTwo = os.path.realpath(sys.argv[summaryTwo])


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

      if "-random1milion" in sys.argv[1:]:
         geckoProfile = 'sourceCode/Python_Scripts/Plot/random_1_milion/random.1milion.profile.xls'
         geckoExonsCount = 'sourceCode/Python_Scripts/Plot/random_1_milion/random.1milion.Exons.Count.txt'
         geckoIntronsCount = 'sourceCode/Python_Scripts/Plot/random_1_milion/random.1milion.Introns.Count.txt'
         geckoPromotersCount = 'sourceCode/Python_Scripts/Plot/random_1_milion/random.1milion.Promoters.Count.txt'
         geckoDNAseCount = 'sourceCode/Python_Scripts/Plot/random_1_milion/random.1milion.DNAse.Count.txt'
         geckoCTCFCount = 'sourceCode/Python_Scripts/Plot/random_1_milion/random.1milion.CTCF.Count.txt'

      subprocess.run(['sourceCode/Python_Scripts/Plot/./radar_chart.py', str(profileFile), str(extProfileFile), str(exonFile), str(intronFile), str(promoterFile), str(dnaseFile),
                     str(ctcfFile), guidesFile, str(mm), str(geckoProfile), str(geckoExonsCount), str(geckoIntronsCount), str(geckoPromotersCount), str(geckoDNAseCount), str(geckoCTCFCount), 
                     str(summaryOne),str(summaryTwo)])

    # -------------------- Error --------------------
    else:
        print("ERROR! \"" + sys.argv[1] + "\" is not an allowed!")

#!/usr/bin/env python

# Python Program to search CRISPR/Cas complex into a genome
import subprocess  # run c++ executable
import time
import os  # instructions manage directories
from os import listdir
from os.path import isfile, join, isdir
import shutil  # remove directory and its content
import sys  # input argv
from subprocess import Popen, PIPE
import glob
import multiprocessing

# path where this file is located
origin_path = os.path.dirname(os.path.realpath(__file__))
# conda path
conda_path = "opt/crispritz/"

VERSION = "2.9.1"


def _run_checked(argv, what):
    """Run a compiled helper and FAIL if it did not exit cleanly.

    subprocess.run() does not raise on a non-zero exit, and nothing here inspected
    the return code, so a crashed binary looked like a completed step. For searchTST
    that is a silent wrong answer rather than a loud one: a segfault writes nothing
    to stderr and leaves a 0-byte targets file, and CRISPRme -- which decides a stage
    failed by looking for output on stderr -- then reports "the search completed
    successfully; the result set is empty" and exits 0. An empty off-target report is
    indistinguishable from a real negative result, so it has to be impossible to
    produce one this way.
    """
    rc = subprocess.run(argv).returncode
    if rc != 0:
        sys.stderr.write(
            f"ERROR: {what} failed ({_exit_detail(rc)}). Its output is incomplete or "
            f"absent, so it must not be treated as an empty result.\n"
        )
        sys.exit(1)


def _exit_detail(rc):
    """'exit code N', or 'signal NAME' for a child killed by a signal."""
    if rc >= 0:
        return f"exit code {rc}"
    import signal as _signal

    try:
        return f"signal {_signal.Signals(-rc).name}"
    except ValueError:
        return f"signal {-rc}"

if "--debug" in sys.argv[1:]:
    # for quick local tests
    corrected_origin_path = origin_path + "/sourceCode/"
else:
    # path corrected to use with conda
    corrected_origin_path = origin_path[:-3] + conda_path


def checkExistance(f_path, element):
    if element == "f":  # check file
        if not isfile(f_path):
            print("ERROR! " + f_path + " is not a file.")
            sys.exit()
    else:  # check directory
        if not isdir(f_path):
            print("ERROR! " + f_path + " is not a directory.")
            sys.exit()


# ALL ACTIVE FUNCTIONS IN CRISPRITZ


def indexGenome():
    if len(sys.argv) < 7 or "help" in sys.argv[1:]:
        print(
            "WARNING: Too few arguments to function index-genome. Please provide:\n",
            "\nEXAMPLE CALL: crispritz.py index-genome name_genome genomeDirectory(FASTA)/ pamFile.txt -bMax 2\n",
            "\n<name_genome>: Name of the genome to create",
            "\n<genomeDirectory>: Directory containing a genome in .fa or .fasta format, need to be separated into single chromosome files.",
            "\n<pamFile>: Text file containing the PAM (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence",
            "\n-bMax <maxBulges_num>: Number of bulges allowed for the search phase",
            "\n-th <num_thread>: (Optional) Number of threads to use. Default uses 4 threads",
        )
        sys.exit()

    nameGenome = sys.argv[2]  # save name of the genome
    dirGenome = os.path.realpath(sys.argv[3])
    checkExistance(dirGenome, "d")
    dirPAM = os.path.realpath(sys.argv[4])
    checkExistance(dirPAM, "f")
    filePAM = open(os.path.realpath(sys.argv[4]), "r")
    listChrs = os.listdir(dirGenome)  # save list of chromosomes

    max_bulges = (sys.argv).index("-bMax") + 1
    try:
        max_bulges = sys.argv[max_bulges]
    except IndexError:
        print("ERROR! No value provided for the -bMax option, please provide a number")
        sys.exit()
    try:
        int(max_bulges)
    except ValueError:
        print(
            "ERROR! The value for the -bMax option is invalid, please provide a number"
        )
        sys.exit()
    # retrive PAM
    PAM = filePAM.read()
    PAM_size = int(PAM.split()[1])
    if PAM_size < 0:
        PAM_size = PAM_size * -1
        PAM = PAM.split()[0][0:PAM_size]
    else:
        PAM = PAM.split()[0][-PAM_size:]

    TSTgenome = (
        PAM + "_" + max_bulges + "_" + nameGenome
    )  # name of the genome in TST format
    dirTSTgenome = "./genome_library/" + TSTgenome  # dir of the genome in TST format

    print(TSTgenome, "Indexing generation:")

    # read number of threads
    th = 4
    if "-th" in sys.argv[1:]:
        try:
            th = (sys.argv).index("-th") + 1
            th = int(sys.argv[th])
            if th < 1:  # if th number is lower than 1 put to 4
                th = 4
        except:
            print(
                "ATTENTION! Check the thread option: -th <th_num> (th_num is an integer)"
            )
            sys.exit()

    if os.path.isdir(dirTSTgenome):  # check if TSTgenome dir exists
        shutil.rmtree(dirTSTgenome)  # remove old TSTgenome dir
    os.makedirs(dirTSTgenome)  # build new TSTgenome dir
    os.chdir(dirTSTgenome)  # move into the TSTgenome dir

    # run buildTST
    start_time = time.time()
    for f in listChrs:
        if f.strip().split(".")[-1] != "fa":
            continue
        print("Indexing:", f)
        _run_checked(
            [
                corrected_origin_path + "buildTST",
                str(dirGenome) + "/" + str(f),
                str(dirPAM),
                str(th),
                max_bulges,
            ],
        "buildTST (index build)",
    )
    print("Finish indexing")
    print("Indexing runtime: %s seconds" % (time.time() - start_time))


def searchTST():
    if len(sys.argv) < 10 or "help" in sys.argv[1:]:
        print(
            "WARNING: Too few arguments to function search with bulges selected\n",
            "\ncrispritz.py search help to see paramaters and example call\n",
        )
        sys.exit()

    nameGenome = os.path.realpath(sys.argv[2])
    PAM = os.path.realpath(sys.argv[3])  # save PAM
    fileGuide = os.path.realpath(sys.argv[4])
    nameResult = sys.argv[5]  # name of result file
    dirTSTgenome = os.path.realpath(sys.argv[2]) + "/"
    try:
        max_bulges = dirTSTgenome.split("/")[-2].split("_")[1]
    except IndexError:
        print(
            "ERROR! The directory name must contain the PAM and the maximum available bulges, please create the Genome Index using crispritz"
        )
        sys.exit()

    checkExistance(PAM, "f")
    checkExistance(fileGuide, "f")
    if not os.path.isdir(dirTSTgenome):  # check if TSTgenome dir exists
        print(
            'ATTENTION! You have to generate the index of "'
            + nameGenome
            + '" with "'
            + PAM
            + '", before the search using index!'
        )
        sys.exit()

    # read number of mismatches
    mm = 0
    if "-mm" in sys.argv[1:]:
        try:
            mm = (sys.argv).index("-mm") + 1
            mm = int(sys.argv[mm])
        except:
            print(
                "ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is a number)"
            )
            sys.exit()

    # read number of bulge RNA
    bRNA = 0
    if "-bRNA" in sys.argv[1:]:
        try:
            bRNA = (sys.argv).index("-bRNA") + 1
            bRNA = int(sys.argv[bRNA])
        except:
            print(
                "ATTENTION! Check the bulge RNA option: -bRNA <bRNA_num> (bRNA_num is a number)"
            )
            sys.exit()

    # read number of bulge DNA
    bDNA = 0
    if "-bDNA" in sys.argv[1:]:
        try:
            bDNA = (sys.argv).index("-bDNA") + 1
            bDNA = int(sys.argv[bDNA])
        except:
            print(
                "ATTENTION! Check the bulge DNA option: -bDNA <bDNA_num> (bDNA_num is a number)"
            )
            sys.exit()

    if int(max_bulges) < bDNA:
        print(
            "WARNING! Max available bulges ("
            + max_bulges
            + ") is smaller than input DNA bulges ("
            + str(bDNA)
            + "). Using "
            + max_bulges
            + " bulges"
        )
        bDNA = max_bulges
    if int(max_bulges) < bRNA:
        print(
            "WARNING! Max available bulges ("
            + max_bulges
            + ") is smaller than input RNA bulges ("
            + str(bRNA)
            + "). Using "
            + max_bulges
            + " bulges"
        )
        bRNA = max_bulges
    # read number of threads
    th = 1
    if "-th" in sys.argv[1:]:
        try:
            th = (sys.argv).index("-th") + 1
            th = int(sys.argv[th])
        except:
            print(
                "ATTENTION! Check the thread option: -th <th_num> (th_num is an integer)"
            )
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

    # Check input correctness
    file_correct_ext = [
        f
        for f in listdir(dirTSTgenome)
        if isfile(join(dirTSTgenome, f)) and not f.endswith(".bin") and not f.startswith(".")
    ]  # Get files not ending with .bin
    if len(file_correct_ext) != 0:  # Some files do not have .bin
        print(
            "ERROR! The directory contains files that are not compatible with the selected search type. Please ensure that the input directory contains only .bin files"
        )
        sys.exit()

    # run searchOnTST
    print("Search START")
    start_time = time.time()
    # optional cap on the TOTAL number of edits (mismatches + DNA + RNA bulges),
    # enforced inside the TST search (CRISPRme issue #107). Absent/negative => disabled.
    max_edits = -1
    if "--max-edits" in sys.argv[1:]:
        try:
            max_edits = int(sys.argv[(sys.argv).index("--max-edits") + 1])
        except:
            print(
                "ATTENTION! Check the total-edits option: --max-edits <n> (n is an integer)"
            )
            sys.exit()
    # optional on-the-fly PAM filter: let ONE generous-flank PAM-less index (an
    # all-N flank of width K, e.g. 8) serve ANY search PAM of length L <= K (e.g.
    # NGG/NRG/NGA/NNGRRT/TTTV) by dropping candidates whose concrete flank does not
    # match the search PAM, before they are emitted. A short search PAM is N-padded
    # up to the stored width K on the distal side. OFF by default; 1 = on. L > K is
    # a hard error (the index never stored those extra flank bases). See README
    # section 3.3.3.
    prune_pam = 1 if "--prune-pam" in sys.argv[1:] else 0
    _run_checked(
        [
            corrected_origin_path + "searchTST",
            str(dirTSTgenome),
            str(fileGuide),
            str(mm),
            str(bDNA),
            str(bRNA),
            str(PAM),
            str(nameResult),
            str(r),
            str(th),
            max_bulges,
            str(max_edits),
            str(prune_pam),
        ],
        "searchTST (index search)",
    )
    print("Search END")
    print("Search runtime: %s seconds" % (time.time() - start_time))


def searchBruteForce():
    if len(sys.argv) < 9 or "help" in sys.argv[1:]:
        print(
            "WARNING: Too few arguments to function search. Please provide:\n",
            "\nEXAMPLE CALL MISMATCHES ONLY: crispritz.py search genomeDirectory/ pamFile.txt guidesFile.txt outputFile -mm 4 -th 4 -t\n",
            "\nEXAMPLE CALL MISMATCHES + BULGES: crispritz.py search indexGenomeDirectory/ pamFile.txt guidesFile.txt outputFile -mm 4 -bDNA 1 -bRNA 2 -th 4 -t\n",
            "\n<genomeDirectory>: Directory containing a genome in .fa or .fasta format (.bin format if bulges present), need to be separated into single chromosome files",
            "\n<pamFile>: Text file containing the PAM sequence (including a number of Ns equal to the guide length) and a space separated number indicating the length of the PAM sequence",
            "\n<guidesFile>: Text file containing one or more guides (including a number of Ns equal to the length of the PAM sequence)",
            "\n<outputFile>: Name of output file",
            "\n-mm <mm_num>: Number of allowed mismatches",
            "\n-bRNA <bRNA_num>: (Optional) Size of RNA bulges",
            "\n-bDNA <bDNA_num>: (Optional) Size of DNA bulges",
            "\n-th < num_thread >: (Optional) Number of threads to use. Default uses all of the available threads (ONE for bulge search)",
            "\n--max-edits <n>: (Optional, indexed search) Cap the TOTAL number of edits (mismatches + DNA bulges + RNA bulges) per off-target. Omitted/negative disables the cap.",
            "\n--prune-pam: (Optional, indexed search) On-the-fly PAM filter. Lets ONE generous-flank PAM-less index (built with an all-N PAM of flank width K, e.g. 8) serve ANY search PAM of length L <= K (e.g. NGG/NRG/NRR/NGA/NNGRRT/TTTV): candidates whose concrete flanking PAM does not match the search PAM are skipped before they are emitted, giving the same off-target loci as a dedicated per-PAM index. A shorter PAM is N-padded up to the stored width K (end for a 3' PAM, front for a 5' PAM). L > K is a hard error. OFF by default (byte-identical to upstream). See README section 3.3.3.",
            "\n{-r,-p,-t}: Output type (-r off-targets list only, -p profile only, -t off-targets AND profile)\n",
        )
        sys.exit()

    genomeDir = os.path.realpath(sys.argv[2]) + "/"
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
                "ATTENTION! Check the mismatches option: -mm <mm_num> (mm_num is an integer)"
            )
            sys.exit()

    # read number of mismatches
    th = 1000
    if "-th" in sys.argv[1:]:
        try:
            th = (sys.argv).index("-th") + 1
            th = int(sys.argv[th])
        except:
            print(
                "ATTENTION! Check the thread option: -th <th_num> (th_num is an integer)"
            )
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

    # Check input correctness
    checkExistance(genomeDir, "d")
    checkExistance(filePAM, "f")
    checkExistance(fileGuide, "f")

    file_correct_ext = [
        f
        for f in listdir(genomeDir)
        if isfile(join(genomeDir, f))
        and not (f.endswith(".fa") or f.endswith(".fasta") or f.endswith(".fai"))
    ]

    if (
        len(file_correct_ext) != 0
    ):  # Some file other than .fa, .fasta or .fai are present in the directory
        print(
            "ERROR! The directory contains files that are not compatible with the selected search type. Please ensure that the input directory contains only .fa files (.fai files are accepted and automatically skipped from the search)"
        )
        sys.exit()

    # run searchBruteForce
    print("Search START")
    start_time = time.time()
    _run_checked(
        [
            corrected_origin_path + "searchBruteForce",
            str(genomeDir),
            str(filePAM),
            str(fileGuide),
            str(mm),
            str(result),
            str(th),
            str(r),
            str(variant),
        ],
        "searchBruteForce (brute-force search)",
    )
    print("Search END")
    print("Search runtime: %s seconds" % (time.time() - start_time))

def _enricher_command():
    """Resolve the enrichment executable, preferring the compiled binary.

    Order: a binary named ``enricher`` on PATH, then the compiled binary shipped
    alongside the reference implementation, else the Python reference script.
    Returns a list suitable to prepend to the 5 positional enricher args.
    """
    found = shutil.which("enricher")
    if found:
        return [found]
    compiled = corrected_origin_path + "Python_Scripts/Enrichment/enricher"
    if os.path.isfile(compiled) and os.access(compiled, os.X_OK):
        return [compiled]
    return [corrected_origin_path + "Python_Scripts/Enrichment/enricher.py"]


def genomeEnrichment_subprocess_VCF(altfile, genfile, dirGenome, doit, dirVCFFiles):
    """Enrich one chromosome; return the enricher's exit status for the parent to
    check. This runs in a Pool worker, so it must not exit the way _run_checked
    does: a SystemExit kills the worker, its task never returns, and pool.join()
    waits forever."""
    return subprocess.run(
        _enricher_command()
        + [
            altfile,
            genfile,
            dirGenome.split("/")[-1],
            str(doit),
            dirVCFFiles,
        ]
    ).returncode


def _enrichment_workers(th_from_flag, th_provided, n_chroms):
    """Choose the number of enrichment workers under a memory budget.

    Each worker loads ~one chromosome, so unbounded fan-out can exhaust RAM.
    Mirrors the CRISPRme post-analysis cap (CRISPRME_MAX_MEM_GB, default 64).
    """
    import os

    try:
        budget = float(os.environ.get("CRISPRME_MAX_MEM_GB", "64"))
    except ValueError:
        budget = 64.0
    try:
        per = float(os.environ.get("CRISPRME_ENRICH_WORKER_GB", "3"))
    except ValueError:
        per = 3.0
    if per <= 0:
        per = 3.0
    mem_cap = max(1, int(budget // per))
    base = th_from_flag if th_provided else (os.cpu_count() or 1)
    return max(1, min(base, mem_cap, n_chroms))


def genomeEnrichment():
    """
    Enrich the genome by parsing the vcf files and then replacing the nucleotides in the reference genome sequence. Compatible with indels
    """
    if len(sys.argv) < 4 or "help" in sys.argv[1:]:
        print(
            "WARNING: Too few arguments to function add-variants. Please provide:",
            "\nEXAMPLE CALL: crispritz.py add-variants vcfFilesDirectory/ genomeDirectory/\n",
            "\n\n<vcfFilesDirectory> : Directory containing VCF files, need to be separated into single chromosome files (multi-sample files will be collapsed into one fake individual)",
            "\n\n<genomeDirectory> : Directory containing a genome in .fa or .fasta format, need to be separated into single chromosome files.",
            "\n\n-th <num_thread>: (Optional) Number of threads to use. Default uses 1 thread",
        )
        sys.exit()

    dirVCFFiles = os.path.realpath(sys.argv[2])
    dirGenome = os.path.realpath(sys.argv[3])
    doit = "no"
    if len(sys.argv) > 4 and sys.argv[4] == "true":
        doit = "yes"
    checkExistance(dirVCFFiles, "d")
    checkExistance(dirGenome, "d")
    # List the per-chromosome VCF files, EXCLUDING tabix/bcftools index sidecars
    # (.tbi/.csi). Build the filtered list in a single comprehension: the previous
    # `for file in listChrs: if file.endswith(".tbi"): listChrs.remove(file)` mutated
    # the list while iterating over it, so it SKIPPED elements and some .tbi files
    # survived. Those survivors were later handed to the enricher (below), which does
    # gzip.open() + readline() on them -- decompressing a binary .tbi index as text
    # raised `UnicodeDecodeError: 'utf-8' codec can't decode byte 0x9d`. See CHANGELOG 2.8.2.
    listChrs = [
        f
        for f in os.listdir(dirVCFFiles)
        if not (f.endswith(".tbi") or f.endswith(".csi"))
    ]

    chr_with_vcf = set()
    # VCF file must contain '.chrN.' with N = number or letter of chr
    for vcf_chr in [
        f
        for f in listdir(dirVCFFiles)
        if isfile(join(dirVCFFiles, f)) and "vcf.gz" in f
    ]:
        for x in vcf_chr.split("."):
            if "chr" in x:
                chr_with_vcf.add(x)
    list_file_ends = [f for f in listdir(dirGenome) if isfile(join(dirGenome, f))]
    file_ends = list_file_ends[0].split(".")[-1]  # For .fa or .fasta
    contains_enr = ""
    if ".enriched." in list_file_ends[0]:
        contains_enr = ".enriched"
    elif ".indels." in list_file_ends[0]:
        contains_enr = ".indels"
    chr_wihtout_vcf = (
        set(
            [
                f.split(contains_enr + ".fa")[0]
                for f in listdir(dirGenome)
                if isfile(join(dirGenome, f))
            ]
        )
        - chr_with_vcf
    )

    # read number of mismatches
    th = 1
    th_provided = "-th" in sys.argv[1:]
    if "-th" in sys.argv[1:]:
        try:
            th = (sys.argv).index("-th") + 1
            th = int(sys.argv[th])
        except:
            print(
                "ATTENTION! Check the thread option: -th <th_num> (th_num must be an integer)"
            )
            sys.exit()

    dirEnrichedGenome = "./variants_genome/"
    if not (os.path.isdir(dirEnrichedGenome)):
        os.makedirs(dirEnrichedGenome)

    os.chdir(dirEnrichedGenome)
    if not (os.path.isdir("./SNPs_genome/")):
        os.makedirs("./SNPs_genome/")
    if not (os.path.isdir("./INDELs_genome/")):
        os.makedirs("./INDELs_genome/")
    os.chdir("./INDELs_genome/")

    memo = open("change_version.txt", "w")

    memo.write(
        "CRISPRitz indels process is now obsolete and has been removed, if you want to process indels you can download our new tool CRISPRme, https://github.com/pinellolab/CRISPRme"
        + "\n"
    )
    memo.write("Thank you")

    memo.close()

    os.chdir("../")
    # os.chdir(dirParsedFiles)
    n_workers = _enrichment_workers(th, th_provided, len(listChrs))
    print(
        "add-variants: enriching %d chromosome(s) with %d worker(s) "
        "(mem budget %s GB, ~%s GB/worker)"
        % (
            len(listChrs),
            n_workers,
            os.environ.get("CRISPRME_MAX_MEM_GB", "64"),
            os.environ.get("CRISPRME_ENRICH_WORKER_GB", "3"),
        ),
        file=sys.stdout,  # benign progress -> stdout, not stderr (stderr trips submit_job [ -s logerror ] gate)
    )
    # Use a 'fork' context: crispritz.py runs its command dispatch at module
    # scope (no __main__ guard), so a 'spawn' start method would re-import and
    # re-invoke add-variants inside every worker. 'fork' is the Linux default
    # (the CRISPRme production platform); we request it explicitly so parallel
    # enrichment behaves identically on macOS.
    try:
        mp_ctx = multiprocessing.get_context("fork")
    except ValueError:
        mp_ctx = multiprocessing
    pool = mp_ctx.Pool(n_workers)

    print("Variants Extraction and Processing START")
    start_time = time.time()
    jobs = []
    for elem in listChrs:
        # os.chdir(dirParsedFiles) #enter in parsed file directory
        chrom = ""
        for cut in elem.split("."):
            if "chr" in cut:
                chrom = cut
        # subprocess.run([corrected_origin_path +
        #                'Python_Scripts/Enrichment/bcf_query.sh', dirVCFFiles+"/"+elem, chrom])
        # altfile = str(chrom + '.alt')
        altfile = dirVCFFiles + "/" + elem
        genfile = str(dirGenome + "/" + chrom + contains_enr + ".fa")
        # pool to process TH vcf file in parallel
        jobs.append((elem, pool.apply_async(
            genomeEnrichment_subprocess_VCF,
            args=(altfile, genfile, dirGenome, doit, dirVCFFiles),
        )))
    # wait until all threads are completed than join
    pool.close()
    pool.join()
    # Nothing used to read these results, so an enricher that crashed, was killed
    # or ran out of disk left add-variants at exit 0 with that chromosome's
    # enriched genome, SNP dict, indel log and contigs missing or truncated -- and
    # CRISPRme's index build, which trusts this exit status, went on to index
    # whatever was there.
    failed = []
    for elem, job in jobs:
        try:
            rc = job.get()
        except Exception as e:  # the worker itself raised
            failed.append((elem, f"{type(e).__name__}: {e}"))
            continue
        if rc != 0:
            failed.append((elem, _exit_detail(rc)))
    if failed:
        for elem, detail in failed:
            sys.stderr.write(
                f"ERROR: enriching {elem} failed ({detail}). Its enriched genome, SNP "
                f"dict, indel log and indel contigs are incomplete or absent.\n"
            )
        sys.exit(1)

    for (
        f
    ) in (
        chr_wihtout_vcf
    ):  # Move chromosomes without vcf to enriched directory and change name adding '.enriched.'
        subprocess.run(
            [
                "cp",
                dirGenome + "/" + f + contains_enr + "." + file_ends,
                "./SNPs_genome/"
                + dirGenome.split("/")[-1]
                + "_enriched/"
                + f
                + ".enriched."
                + file_ends,
            ]
        )

    print("Variants Extraction and Processing END")
    print("Runtime: %s seconds" % (time.time() - start_time))


def removeFile(to_remove):
    """
    Remove file given his path in input
    """
    try:
        os.remove(to_remove)
    except:
        pass


# return version


def version():
    print("CRISPRitz v", VERSION)


# HELP FUNCTION
def callHelp():
    print("version", VERSION)
    print(
        "help:\n",
        "\nALL FASTA FILEs USED BY THE SOFTWARE MUST BE UNZIPPED AND CHROMOSOME SEPARATED, ALL VCFs USED BY THE SOFTWARE MUST BE ZIPPED AND CHROMOSOME SEPARATED",
        "\n",
        "\ncrispritz.py add-variants FUNCTION TO ADD VARIANTS DATA TO A FASTA GENOME",
        "\ncrispritz.py index-genome FUNCTION TO CREATE GENOME INDEX TO PERFORM FAST SEARCHES WITH BULGES",
        "\ncrispritz.py search FUNCTION TO PERFORM SEARCHES ON A GENOME (INDEXED OR PLAIN FASTA)",
        "\n",
        "\nADD help TO ANY FUNCTION TO VISUALIZE A BRIEF HELP PAGE (example: crispritz.py index-genome help)\n",
    )


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
elif sys.argv[1] == "add-variants":
    genomeEnrichment()
else:
    print('ERROR! "' + sys.argv[1] + '" is not an allowed!')

"""
"""

from utils import INDELS_GENOME, SNP_GENOME, VARIANTS_GENOME

from time import time

import multiprocessing
import signal
import sys
import os


# TODO: use argparse for args
def enrich_genome(vcf: str, genome: str, cores: int, outdir: str) -> None:
    """
    """

    isvcfdir = True if os.path.isdir(vcf) else False  # directory containing VCFs or single VCF
    isgenomedir = True if os.path.isdir(genome) else False  # directory containing genome or single file
    assert os.path.exists(vcf)
    assert os.path.exists(genome)
    assert os.path.exists(outdir)
    if isvcfdir:  # list all the vcfs in the directory
        vcf = [f for f in os.listdir(vcf) if ".vcf" in f and not f.endswith(".tbi")]
    # create the enriched genome directory tree
    enriched_genome_dir = os.path.join(outdir, VARIANTS_GENOME)
    if not os.path.isdir(enriched_genome_dir):
        os.mkdir(enriched_genome_dir)
    snps_genome_dir = os.path.join(enriched_genome_dir, SNP_GENOME)
    if not os.path.isdir(snps_genome_dir):
        os.mkdir(snps_genome_dir)
    indels_genome_dir = os.path.join(enrich_genome, INDELS_GENOME)
    if not os.path.isdir(indels_genome_dir):
        os.mkdir(indels_genome_dir)
    # capture SIGINT with multiprocessing Pool
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    # create threads pool
    pool = multiprocessing.Pool(processes=cores)
    signal.signal(signal.SIGINT, original_sigint_handler)
    sys.stderr.write("Variants extraction and processing START\n")
    start = time()
    stop = time()
    sys.stderr.write("Variants Extraction and Processing END")
    sys.stderr.write("Runtime: %s seconds\n" % ((stop - start)))
    






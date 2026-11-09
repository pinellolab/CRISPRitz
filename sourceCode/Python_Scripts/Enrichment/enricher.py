#!/usr/bin/env python

import sys
import re
import time
import textwrap
import os
import json
import pandas as pd
import gzip

start_time = time.time()

print("READING VCF FILE AND CHROMOSOME")

altFile = sys.argv[1]  # file with variations
genomeFile = sys.argv[2]  # genome file
dir_enr_name = sys.argv[3] + '_enriched' #name of directory for saving files
inAltFile = gzip.open(altFile, "rt")  # variations file open
inGenomeFile = open(genomeFile, "r")  # genome file open
vcfName = os.path.basename(sys.argv[5])

#read fasta header (chr_name)
genomeHeader = inGenomeFile.readline()
currentChr = genomeHeader[1:].strip()
if 'chr' not in currentChr:
    currentChr = 'chr'+currentChr

#read VCF header
VCFheader=inAltFile.readline().split('\t')

#all the fasta content is read and stored in a single line with \n removed and uppercased to avoid possible conflicts
genomeStr = inGenomeFile.read()
genomeStr = genomeStr.replace('\n', '')
genomeStr = genomeStr.upper()

# Q is the symbol written where a position's SNP ALTs, with its reference base when that
# is A, C, G or T, cover all four bases (SNPsProcess), and buildTST, searchTST and
# searchBruteForce read it as any base: a FASTA that already holds one would pass a base
# of its own off as such a position.
if 'Q' in genomeStr:
    sys.stderr.write("ERROR: %s holds 'Q' at %d: a reference FASTA must not hold Q, which is "
                     "reserved for four-allele SNP positions\n"
                     % (currentChr, genomeStr.index('Q') + 1))
    sys.exit(1)

#genome str is listed to be modified faster and easier
genomeList = list(genomeStr)

# The A/C/G/T bases coded so far at each position, as a 4-bit mask (0 = no SNP record
# here yet). A multi-allelic site normalised with `bcftools norm -m -any` arrives as one
# record per ALT, so each record must add its ALT to what the earlier ones coded instead
# of rebuilding the code from the reference base. The mask is kept beside genomeList
# because a letter there cannot tell a code this script wrote from an ambiguity code the
# reference FASTA already holds.
snp_mask = bytearray(len(genomeStr))
# positions coded Q (all four bases), counted for the stdout note
four_base_count = 0

iupac_code = {}

iupac_code.setdefault('R', []).append('AG')
iupac_code.setdefault('R', []).append('GA')

iupac_code.setdefault('Y', []).append('CT')
iupac_code.setdefault('Y', []).append('TC')

iupac_code.setdefault('S', []).append('GC')
iupac_code.setdefault('S', []).append('CG')

iupac_code.setdefault('W', []).append('AT')
iupac_code.setdefault('W', []).append('TA')

iupac_code.setdefault('K', []).append('GT')
iupac_code.setdefault('K', []).append('TG')

iupac_code.setdefault('M', []).append('AC')
iupac_code.setdefault('M', []).append('CA')

iupac_code.setdefault('B', []).append('CGT')
iupac_code.setdefault('B', []).append('GCT')
iupac_code.setdefault('B', []).append('TGC')
iupac_code.setdefault('B', []).append('GTC')
iupac_code.setdefault('B', []).append('CTG')
iupac_code.setdefault('B', []).append('TCG')

iupac_code.setdefault('D', []).append('AGT')
iupac_code.setdefault('D', []).append('GAT')
iupac_code.setdefault('D', []).append('TAG')
iupac_code.setdefault('D', []).append('ATG')
iupac_code.setdefault('D', []).append('GTA')
iupac_code.setdefault('D', []).append('TGA')

iupac_code.setdefault('H', []).append('ACT')
iupac_code.setdefault('H', []).append('CAT')
iupac_code.setdefault('H', []).append('TCA')
iupac_code.setdefault('H', []).append('ATC')
iupac_code.setdefault('H', []).append('CTA')
iupac_code.setdefault('H', []).append('TAC')

iupac_code.setdefault('V', []).append('ACG')
iupac_code.setdefault('V', []).append('CAG')
iupac_code.setdefault('V', []).append('GAC')
iupac_code.setdefault('V', []).append('AGC')
iupac_code.setdefault('V', []).append('CGA')
iupac_code.setdefault('V', []).append('GCA')

iupac_code.setdefault('N', []).append('ACGT')
iupac_code.setdefault('N', []).append('CAGT')
iupac_code.setdefault('N', []).append('GACT')
iupac_code.setdefault('N', []).append('AGCT')
iupac_code.setdefault('N', []).append('CGAT')
iupac_code.setdefault('N', []).append('GCAT')
iupac_code.setdefault('N', []).append('GCTA')
iupac_code.setdefault('N', []).append('CGTA')
iupac_code.setdefault('N', []).append('TGCA')
iupac_code.setdefault('N', []).append('GTCA')
iupac_code.setdefault('N', []).append('CTGA')
iupac_code.setdefault('N', []).append('TCGA')
iupac_code.setdefault('N', []).append('TAGC')
iupac_code.setdefault('N', []).append('ATGC')
iupac_code.setdefault('N', []).append('GTAC')
iupac_code.setdefault('N', []).append('TGAC')
iupac_code.setdefault('N', []).append('AGTC')
iupac_code.setdefault('N', []).append('GATC')
iupac_code.setdefault('N', []).append('CATG')
iupac_code.setdefault('N', []).append('ACTG')
iupac_code.setdefault('N', []).append('TCAG')
iupac_code.setdefault('N', []).append('CTAG')
iupac_code.setdefault('N', []).append('ATCG')
iupac_code.setdefault('N', []).append('TACG')

iupac_code_scomposition = {
          "R":"AG",
          "Y":"CT",
          "S":"GC",
          "W":"AT",
          "K":"GT",
          "M":"AC",
          "B":"CGT",
          "D":"AGT",
          "H":"ACT",
          "V":"ACG",
          "r":"AG",
          "y":"CT",
          "s":"GC",
          "w":"AT",
          "k":"GT",
          "m":"AC",
          "b":"CGT",
          "d":"AGT",
          "h":"ACT",
          "v":"ACG",
          "A":"A",
          "T":"T",
          "C":"C",
          "G":"G",
          "a":"A",
          "t":"T",
          "c":"C",
          "g":"G",
          'N':'ATGC'
        }


# base SET -> IUPAC letter, order-independent: the A/C/G/T bases present as a 4-bit
# mask. Built from iupac_code, which has no single-base key, so one base maps to
# nothing and leaves the genome unchanged. All four bases are Q, not N: N is what
# buildTST and searchTST read as an assembly gap, which matches nothing, so writing it
# removed every target through the position from the index, the reference one too.
_BASE_BIT = {'A': 1, 'C': 2, 'G': 4, 'T': 8}
_IUPAC_BY_MASK = {sum(_BASE_BIT[c] for c in spellings[0]): letter
                  for letter, spellings in iupac_code.items()}
_IUPAC_BY_MASK[15] = 'Q'
# the reference letters a SNP is coded over: A/C/G/T and the IUPAC ambiguity codes. N is
# a gap, and so is any other byte (buildTST and searchBruteForce read it as N); Q is
# refused when the FASTA is read.
_REFERENCE_CODES = frozenset('ACGTRYSWKMBDHV')


# A plain decimal number (what an AF value is). Spelled out rather than left to float(),
# which also takes 'nan', 'inf', '1_000' and surrounding spaces, so that enricher.cpp
# (strtod behind the same check) reads every AF string the same way.
_AF_NUMBER = re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z')


def _af_rank(af):
    """An AF value string as a number, to tell which of two records of one allele is the
    more frequent (_dict_put): 0 when it is missing or not a plain decimal number."""
    return float(af) if _AF_NUMBER.match(af) else 0.0


def SNPsProcess(line):
    global four_base_count
    if len(line[3]) > 1: #se len ref > 1 significa che non è snp
        return
    line[1] = str(int(line[1])-1)  # taaac per corregere conta posizione da zero/uno
    pos = int(line[1])
    if not 0 <= pos < len(genomeStr):
        # A POS past the end raised IndexError here; POS 0 read and coded the LAST base
        # (index -1), and enricher.cpp read and wrote out of bounds. Fail, both ports alike.
        sys.stderr.write("ERROR: %s:%d is outside the chromosome (%d bp): is the VCF on "
                         "another assembly?\n" % (currentChr, pos + 1, len(genomeStr)))
        sys.exit(1)
    altAlleles = line[4].strip().split(',') #list di tutti gli alleli alternative
    referenceNucleotide = genomeStr[pos] #nucleotide reference preso dalla stringa del fasta
    if referenceNucleotide not in _REFERENCE_CODES:
        # An N is an assembly gap, which buildTST and searchTST match to nothing, and a
        # variant record does not make it a base: the enriched genome holds exactly the
        # reference's N (2.9.1 coded a record listing two or three ALTs over an N, making
        # a gap searchable, and merging split records would have too). Any other byte
        # that is no base or IUPAC code (such as the X of RepeatMasker -x) is a gap too,
        # which buildTST and searchBruteForce read as N, and stays as it is. The dict
        # still lists the record.
        return
    # The reference base and every single-character ALT, as A/C/G/T bits. A member
    # that is not A/C/G/T -- the '*' of a spanning deletion, an 'N', a lowercase
    # base, an IUPAC reference base -- adds nothing, so the real SNP alleles beside
    # it are still coded. Building the code from the whole set instead left the base
    # unchanged whenever such a member was present, dropping a real SNP at every
    # multiallelic '*' site. This is enricher.cpp's rule (baseBit/iupacFromMask).
    #
    # The set starts from what earlier records coded HERE, not from the reference base
    # alone. A normalised VCF (`bcftools norm -m -any`, how every CRISPRme panel is
    # built) splits a multi-allelic site into one record per ALT, and rebuilding the
    # code from the reference on each record kept only the LAST record's ALT in the
    # file: the others reached neither the enriched genome nor the search. Now the
    # code is that of the site written as one record with every ALT, in any order;
    # all four bases are Q (_IUPAC_BY_MASK), which the search reads as any base.
    mask = snp_mask[pos] or _BASE_BIT.get(referenceNucleotide, 0)
    for alt in altAlleles:
        if len(alt) == 1:
            mask |= _BASE_BIT.get(alt, 0)
    if mask == 15 and snp_mask[pos] != 15:
        four_base_count += 1
    snp_mask[pos] = mask
    letter = _IUPAC_BY_MASK.get(mask)
    if letter is not None:
        genomeList[pos] = letter #inserisco il IUPAC nt nel genoma enriched


# def SNPsProcess(line):
#     #x = line.strip().split('\t') #split alt file to add snps to genome
#     x = line
#     #del x[1] #remove the unnecessary rsID from the split list
#     x[1] = str(int(x[1])-1)  # taaac
#     if (',' in x[4]) and (len(x[3]) == 1) and ('>' not in x[4]) and (len(x[4]) < 6):
#         altstr = str(x[4])
#         k = altstr.split(',')
#         if (len(x[4]) == 3):
#             if (str(genomeList[int(x[1])])) not in k:
#                 original = genomeList[int(x[1])]
#             else:
#                 original = ''
#             snp = k[0]+k[1]
#             iupacvalue = str(original+snp)  #TODO if original is a single nucleotide, use iupac_code_scomposition[original]
#             for key, value in iupac_code.items():
#                 if iupacvalue in value:
#                     genomeList[int(x[1])] = str(key)
#         elif (len(x[4]) == 5):
#             if (len(k[0]) == 1) and (len(k[1]) == 1) and (len(k[2]) == 1):
#                 if (str(genomeList[int(x[1])])) not in k:
#                     original = genomeList[int(x[1])]
#                 else:
#                     original = ''
#                 original = genomeList[int(x[1])]
#                 snp = k[0]+k[1]+k[2]
#                 iupacvalue = str(original+snp)  #TODO if original is a single nucleotide, use iupac_code_scomposition[original]
#                 for key, value in iupac_code.items():
#                     if iupacvalue in value:
#                         genomeList[int(x[1])] = str(key)
#     elif (',' not in x[4]) and (len(x[3]) == 1) and ('>' not in x[4]) and (len(x[4]) == 1):
#         original = iupac_code_scomposition[genomeList[int(x[1])]]
#         snp = x[4]
#         iupacvalue = str(original+snp)
#         for key, value in iupac_code.items():
#             if iupacvalue in value:
#                 genomeList[int(x[1])] = str(key)

def chromosomeSave():
    genomeStr = "".join(genomeList)
    os.chdir("./SNPs_genome/")
    # add-variants runs one enricher per chromosome in parallel, all into this one
    # directory: two that reach here together both saw it missing and the second
    # mkdir raised, losing that chromosome's outputs. enricher.cpp treats EEXIST as OK.
    os.makedirs(dir_enr_name, exist_ok=True)
    outFile = open(dir_enr_name + '/' + genomeHeader[1:(len(genomeHeader)-1)]+'.enriched'+'.fa', 'w')
    outFile.write(genomeHeader+genomeStr+'\n')
    outFile.close()

def _gt_alt_alleles(gt):
    """ALT allele indices carried by a VCF GT field.

    Splits on '|' and '/', so phased, unphased, haploid ('1') and
    partially-missing ('.|1') genotypes are all handled. Reference (0) and
    missing ('.') alleles are skipped, so the result holds only ALT indices.
    Matching on these integers -- instead of searching the raw GT text for an
    allele number -- is what keeps '0|12' from matching allele 1.
    """
    alleles = set()
    for a in gt.replace('/', '|').split('|'):
        if a.isdigit():
            v = int(a)
            if v > 0:  # 0 is the reference allele
                alleles.add(v)
    return alleles


def _af_value(line):
    # Raw INFO/AF value string of this record, found by its key, or "" if the record
    # has none. "" is the deliberate placeholder (not "."): it's the only value besides
    # "0" that pandas.to_numeric() parses as NaN instead of raising downstream. The key
    # must be exactly "AF" (not a prefix match, which would take a gnomAD-style AF_afr).
    # Until 2.9.1 the AF field's place was taken from the FIRST record only, so every
    # record whose INFO lists its fields in another order or number got "": a merged VCF
    # lost the AF of 224,981 of 232,384 chr22 indels that way, and HPRC chrY 91,950 of
    # its 102,665 indel AFs.
    for entry in line[7].split(";"):
        key, sep, value = entry.partition("=")
        if key == "AF":
            return value
    return ""


def _entry_refalt(entry):
    """The "<REF>,<ALT>" of one dict entry "<samples>;<REF>,<ALT>;<rsID>;<AF>", read from
    the right, so a ';' in a sample name cannot move it (the rsID's ';' is already ',')."""
    return entry.rsplit(';', 3)[1]


def _gt_for_allele(gt, k):
    """GT ``gt`` with allele k written 1 and every other called allele 0 (the REF or
    another ALT), missing ones and separators as they were: the form a split record's
    genotypes have, where 1 is its one ALT."""
    return ''.join(('1' if int(a) == k else '0') if a.isdigit() else a
                   for a in re.split(r'([|/])', gt))


def _gt_union(gt1, gt2):
    """One sample's genotypes of one allele on two records, the earlier first, both in
    split form (1 = the allele). Phased alike (as many alleles, every separator '|'): a
    slot holds the allele if either record puts it there, else a called 0 if either has
    one. Otherwise, as unphased slots cannot be matched up, the genotype with more copies
    of the allele (on a tie, the earlier)."""
    t1, t2 = re.split(r'([|/])', gt1), re.split(r'([|/])', gt2)
    if len(t1) != len(t2) or '/' in t1 or '/' in t2:
        return gt2 if t2.count('1') > t1.count('1') else gt1
    out = []
    for x, y in zip(t1, t2):
        if x == '|':
            out.append(x)
        elif '1' in (x, y):
            out.append('1')
        elif '0' in (x, y):
            out.append('0')
        else:
            out.append(x)
    return ''.join(out)


def _merge_entries(old, k_old, new, k_new):
    """One allele's dict entries from two records (the same REF>ALT twice): every carrier
    of either, with genotypes in split form (_gt_for_allele) and a carrier of both
    getting the union of its two (_gt_union); the rsID and AF of the more frequent
    record (highest INFO/AF; on a tie, the earlier). Until 2.9.1 the later record's
    entry replaced the earlier, dropping the earlier record's carriers."""
    s_old, refalt, rsid, af = old.rsplit(';', 3)
    s_new, _, rsid_new, af_new = new.rsplit(';', 3)
    gts = {}
    for samples, k in ((s_old, k_old), (s_new, k_new)):
        for carrier in samples.split(',') if samples else []:
            name, gt = carrier.rsplit(':', 1)
            gt = _gt_for_allele(gt, k)
            gts[name] = _gt_union(gts[name], gt) if name in gts else gt
    if _af_rank(af_new) > _af_rank(af):
        rsid, af = rsid_new, af_new
    return ','.join(sorted(n + ':' + g for n, g in gts.items())) + ';' + refalt + ';' + rsid + ';' + af


def _dict_put(key, entries, alleles):
    """Store one record's per-ALT entries under "chr,pos", MERGED with what earlier
    records stored there; alleles[i] is the allele index of entries[i]'s ALT in its
    record (1 on a split record).

    A position's value is its ALTs' entries joined by '$', the form a record listing
    several ALTs has always produced, so every reader of multi-ALT values reads it. Each
    entry keeps its own record's carriers, genotypes, rsID and AF: the genotypes of a
    split record use 1 for its ALT, as tier0_registry's default alt_index expects.
    Assigning instead let a later record at the position overwrite the earlier ones, so
    a split multi-allelic site kept only its last ALT. An allele already listed here
    (the same REF>ALT on two records) gets the carriers of both (_merge_entries); any
    other ALT is appended, in file order."""
    old = chr_dict_snps.get(key)
    if old is None:
        chr_dict_snps[key] = '$'.join(entries)
        for e, k in zip(entries, alleles):
            if k != 1:
                entry_allele[(key, _entry_refalt(e))] = k
        return
    merged = old.split('$')
    refalts = [_entry_refalt(e) for e in merged]
    for e, k in zip(entries, alleles):
        ra = _entry_refalt(e)
        if ra in refalts:
            j = refalts.index(ra)
            merged[j] = _merge_entries(merged[j], entry_allele.pop((key, ra), 1), e, k)
        else:
            merged.append(e)
            refalts.append(ra)
            if k != 1:
                entry_allele[(key, ra)] = k
    chr_dict_snps[key] = '$'.join(merged)


def add_to_dict_snps(line):
    list_samples = []
    list_chars = []
    if len(line[3]) == 1 and len(line[4]) == 1:
        for pos, i in enumerate(line[9:]):          #if sample has 1|1 0|1 or 1|0, #NOTE may change for different vcf
            if ('1' in i.split(':')[0]):
                list_samples.append(VCFheader[ pos + 9]+':'+i.split(':')[0])
        
        # print(list_samples)
        
        chr_pos_string = currentChr + ',' + line[1] #chr,position
        #Add in last two position the ref and alt nucleotide, eg: chrX,100 -> sample1,sample5,sample10;A,T;rsID100;0.01
        #If no sample was found, the dict is chrX,100 -> ;A,T;rsID100;0.01
        rsID = line[2]
        list_chars.append(line[3])
        list_chars.append(line[4])
        af = _af_value(line)
        if len(list_samples) > 0:
            _dict_put(chr_pos_string, [','.join(sorted(list_samples)) + ';' + ','.join(list_chars) + ";" + rsID + ";" + af], [1])
        else:
            _dict_put(chr_pos_string, [';' + ','.join(list_chars) + ";" + rsID + ";" + af], [1]) #None
        # print(chr_dict_snps[chr_pos_string])
    elif len(line[3]) == 1:
        variants = line[4].split(",")
        snps = []
        values_for_allele_info = []
        for value, var in enumerate(variants):
            if len(var) == 1: 
                snps.append(var)
                values_for_allele_info.append(value+1) #save number corresponding to snp (1,2,3...)
        dict_of_lists_samples = {}
        for snp in snps:
            dict_of_lists_samples[snp] = []
        if len(snps) > 0:		
            for pos, sample in enumerate(line[9:]):
                alleles = _gt_alt_alleles(sample.split(':')[0])
                # no break: a sample carrying two different ALTs belongs to BOTH
                for idx, value in enumerate(values_for_allele_info):
                    if value in alleles:
                        dict_of_lists_samples[snps[idx]].append(VCFheader[ pos + 9]+':'+sample.split(':')[0]) #add to correct entry of dict the sample with such snp
                    
            chr_pos_string = currentChr + ',' + line[1]
            rsID = line[2].split(',')
            af = _af_value(line).split(',')

            final_entry = []
            for idx, snp in enumerate(snps):
                list_chars = [line[3]]
                list_chars.append(snp)
                af_idx = values_for_allele_info[idx] - 1
                af_val = af[af_idx] if af_idx < len(af) else ""  # multiallelic site may have fewer AF values than ALT alleles
                if len(dict_of_lists_samples[snp]) > 0:
                    final_entry.append(','.join(sorted(dict_of_lists_samples[snp])) + ';' + ','.join(list_chars) + ";" + rsID[0] + ";" + af_val)
                else:
                    final_entry.append(';' + ','.join(list_chars) + ";" + rsID[0] + ";" + af_val)
            _dict_put(chr_pos_string, final_entry, values_for_allele_info)

def dictSave():
    #os.chdir("./SNPs_genome/")
    with open('my_dict_' + currentChr + '.json', 'w') as f:
        json.dump(chr_dict_snps, f) 


# The fake-indel window (enricher.cpp's INDEL_LEFT / INDEL_RIGHT): each indel allele's
# contig is INDEL_LEFT reference bases, then the record's REF (from POS on) with the ALT
# written in its place, then INDEL_RIGHT reference bases (fewer where the chromosome
# ends); one N stands between two contigs. A target of L bases (guide + PAM) with b DNA
# bulges holds a record's change only where L + b <= INDEL_LEFT + 2 and
# L + b <= INDEL_RIGHT + 1 (INDEL_LEFT + 1 for a record whose REF and ALT differ at their
# first base), so LEFT = RIGHT = Lmax - 1 + bmax: Lmax = 29 (a 21-nt guide and a PAM of 8
# nt, the longest the general K=8 index serves), bmax = 2 (the most DNA bulges CRISPRme's
# web offers). Up to 2.9.1 the
# window was 25 / 27 (L + b <= 27). The log's CHR field records each window's real start,
# from which a reader derives the left flank (POS - 1 - start).
INDEL_LEFT = 30
INDEL_RIGHT = 30


def indel_to_fasta(line, id_indel, start_fake_pos):
    if (len(line[3]) > 1 or len(line[4]) > 1) and '<' not in line[3]:
        #print(''.join(line[3:5]))
        indels = []
        values_for_allele_info = []
        if ',' in line[4]:
            splitted = line[4].split(',')
            for value, s in enumerate(splitted):
                if (len(s) != len(line[3]) or (len(s) > 1 and len(line[3]) > 1)) and '<' not in s:
                    indels.append(s)
                    values_for_allele_info.append(value+1)
        elif (len(line[4]) != len(line[3]) or (len(line[4]) > 1 and len(line[3]) > 1)) and '<' not in line[4]:
            indels.append(line[4])
            values_for_allele_info.append(1)
        #print(search_sample_value, line[3], line[4])
        if len(indels) > 0:

            dict_of_lists_samples = {}
            for indel in indels:
                dict_of_lists_samples[indel] = []

            for pos, sample in enumerate(line[9:]):          #if sample has 1|1 0|1 or 1|0, #NOTE may change for different vcf
                alleles = _gt_alt_alleles(sample.split(':')[0])
                # no break: a sample carrying two different ALTs belongs to BOTH
                for idx, value in enumerate(values_for_allele_info):
                    if value in alleles:
                        dict_of_lists_samples[indels[idx]].append(VCFheader[ pos + 9]) #add to correct entry of dict the sample with such snp
            if len(indels) > 0:
                rsID = line[2].split(',')
                af = _af_value(line).split(',')

                start_position = int(line[1]) - 1 - INDEL_LEFT
                end_position = int(line[1]) - 1 + len(line[3]) + INDEL_RIGHT
                if start_position < 0:
                    # The window is anchored so the REF base lands at index INDEL_LEFT,
                    # which a negative start cannot satisfy: Python would slice from the
                    # far end of the chromosome and hand back an EMPTY window, so the
                    # contig and the coordinates logged for it would disagree.
                    # (enricher.cpp throws instead, so the two ports only agree if both
                    # skip.)
                    print(f"WARNING: skipping indel at {currentChr}:{line[1]} -- fewer than {INDEL_LEFT + 1} bases from the chromosome start")
                    return id_indel, start_fake_pos
                sub_fasta = genomeStr[start_position:end_position]
                for idx, indel in enumerate(indels):
                    if len(dict_of_lists_samples[indel]) > 0:
                        indel_info = f"{currentChr}_{line[1]}_{line[3]}_{indel}"
                        
                        #sub_fasta[25] = line[3] 
                        #sub_fasta = ''.join(sub_fasta)
                        # Build THIS allele's window from the reference window every
                        # time. Assigning back into sub_fasta fed the previous allele's
                        # mutated sequence into this substitution, so a site carrying two
                        # ALTs emitted a contig holding BOTH of them -- a haplotype no
                        # sample need carry -- and its length, refseq and FAKEPOS drifted
                        # with it. Split (biallelic) input hides this, because the loop
                        # runs once per record; input merged into multiallelic sites (see
                        # gnomAD-converter --multiallelic) does not. Dropping the `break`
                        # in the carrier scan made it worse, by letting more alleles reach
                        # this gate.
                        allele_fasta = sub_fasta[0:INDEL_LEFT] + re.sub(line[3], indel, sub_fasta[INDEL_LEFT:], 1, flags=re.IGNORECASE)

                        #fasta_out.write(f'>{currentChr}_{start_position}-{end_position}_{id_indel}\n')
                        list_fasta_indels.append(allele_fasta+'\n'+"N\n")
                        
                        refseq = genomeStr[start_position:start_position+len(allele_fasta)]
                        end_fake_pos = start_fake_pos + len(allele_fasta)#(end_position - start_position)

                        af_idx = values_for_allele_info[idx] - 1
                        af_val = af[af_idx] if af_idx < len(af) else ""  # multiallelic site may have fewer AF values than ALT alleles -- "" (not ".") so pd.to_numeric() reads it as NaN instead of crashing
                        log_indels.append([f"{currentChr}_{start_position}-{end_position}_{id_indel}", ",".join(dict_of_lists_samples[indel]), rsID[0], af_val, indel_info, f"{start_fake_pos},{end_fake_pos}", refseq])
                    
                        id_indel += 1
                        start_fake_pos = end_fake_pos + 1 

    return id_indel, start_fake_pos

def logIndelsSave():
    fasta_out.write(''.join(list_fasta_indels))
    fasta_out.close()
    df = pd.DataFrame(log_indels, columns=['CHR','SAMPLES','rsID','AF','indel', 'FAKEPOS', 'refseq'])
    df.to_csv('log' + currentChr + '.txt', index=False, sep='\t')
    #dict_indels = df.to_dict(orient='index')
    #with open('log' + currentChr + '.json', 'w') as f:
    #    json.dump(log_indels, f) 

print('START ENRICHMENT WITH SNVs AND SVs')

if sys.argv[4] == 'yes':
    chr_dict_snps = dict()
    # (key, "REF,ALT") -> the allele index of that entry's ALT in its record, for an
    # entry of a record listing several ALTs (1, a split record's, is not stored)
    entry_allele = dict()
    log_indels = []
    id_indel = 1
    start_fake_pos = 0
    if not os.path.isdir(f"fake_{vcfName}_{currentChr}"):
        os.mkdir(f"fake_{vcfName}_{currentChr}")
    output_dir_indels = os.path.abspath(f"fake_{vcfName}_{currentChr}")
    fasta_out = open(f'{output_dir_indels}/fake{currentChr}.fa', 'w+')
    fasta_out.write(f">fake{currentChr}\n")
    list_fasta_indels = []

for line in inAltFile:
    if ('#CHROM') in line:
        VCFheader = line.strip().split('\t')   #Save this header for retrieving sample id
        break
for line in inAltFile:
    line = line.strip().split('\t')
    # issue #143: the variant-info string embedded per position is ';'-delimited
    # (samples;ref,alt;rsID;AF). The VCF ID column (line[2]) can itself contain
    # ';' for dbSNP multi-rsID records (e.g. "rs1;rs2"), which collides with that
    # delimiter and shifts AF (and every field after) by one -> downstream
    # float(rsID) crashes. Normalize the ID field's ';' to ',' (CRISPRme's
    # multi-value convention) once, so all rsID uses below stay a single field.
    if len(line) > 2:
        line[2] = line[2].replace(";", ",")
    if line[6] not in ('PASS', '.'):  # '.' = no filter applied (e.g. HPRC vcfwave VCFs)
        continue
    if sys.argv[4] == 'yes': #if true do all the creation, if false do only genome enrichment with SNPs
        #print(sys.argv[4])
        add_to_dict_snps(line)
        id_indel, start_fake_pos = indel_to_fasta(line, id_indel, start_fake_pos)
    SNPsProcess(line)


if four_base_count:
    print("NOTE: %d position(s) of %s are coded Q, the four-allele symbol, which the search "
          "reads as any base (N stays a gap)" % (four_base_count, currentChr))
#saving chr after enrichment with snps
chromosomeSave()
if sys.argv[4] == 'yes':
    dictSave()
    logIndelsSave()

#REMOVED INDELS CREATION AND SUBSTITUTED WITH NEW SCRIPT

#step to generate fasta files with indels and snps, now INDELS are treated separately, this works only for single sample VCF
# for line in inAltFile:
#     # if 'PASS' not in line: #skip line with QUAL inferior to PASS
#     #     continue
#     x = line.strip().split('\t')
#     x[0] = str(int(x[0])-1)  # taaac
#     if (',' not in x[2]) and (',' not in x[1]) and ('>' not in x[2]) and (len(x[1]) == 1) and (len(x[2]) > 1):
#         genomeList[int(x[0])] = str(x[2])
#     elif (',' not in x[2]) and (',' not in x[1]) and ('>' not in x[2]) and (len(x[1]) > 1) and (len(x[2]) == 1):
#         point = int(x[0])
#         fine = point+len(x[1])
#         genomeList[point] = str(x[2])
#         genomeList[point+1:fine] = [''] * (len(x[1])-1)

# genomeStr = "".join(genomeList)

# os.chdir("../")
# os.chdir("./INDELs_genome/")
# if not (os.path.isdir(dir_enr_name)):
#     os.mkdir(dir_enr_name)
# outfile = open(dir_enr_name + '/' + genomeHeader[1:(len(genomeHeader)-1)]+'.indels'+'.fa', 'w')
# outfile.write(genomeHeader+genomeStr+'\n')

print("DONE IN ",int(time.time()-start_time),' seconds')

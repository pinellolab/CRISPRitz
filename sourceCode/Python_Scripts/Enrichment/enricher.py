#!/usr/bin/env python

import sys
import re
import time
import textwrap
import os
import json
import pandas as pd
import gzip
from array import array

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

#genome str is listed to be modified faster and easier
genomeList = list(genomeStr)

# The A/C/G/T bases coded so far at each position, as a 4-bit mask (0 = no SNP record
# here yet). A multi-allelic site normalised with `bcftools norm -m -any` arrives as one
# record per ALT, so each record must add its ALT to what the earlier ones coded instead
# of rebuilding the code from the reference base. The mask is kept beside genomeList
# because a letter there cannot tell a code this script wrote from an ambiguity code the
# reference FASTA already holds.
snp_mask = bytearray(len(genomeStr))
# Every single-base A/C/G/T ALT seen, with its frequency (_rank_af), for the positions
# that end up holding all four bases (see resolveFourBasePositions). ~17 bytes per ALT.
alt_af_pos = array('q')
alt_af_bit = array('B')
alt_af_val = array('d')
four_base_positions = set()

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
# nothing and leaves the genome unchanged.
_BASE_BIT = {'A': 1, 'C': 2, 'G': 4, 'T': 8}
_IUPAC_BY_MASK = {sum(_BASE_BIT[c] for c in spellings[0]): letter
                  for letter, spellings in iupac_code.items()}


# A plain decimal number (what an AF value is). Spelled out rather than left to float(),
# which also takes 'nan', 'inf', '1_000' and surrounding spaces, so that enricher.cpp
# (strtod behind the same check) reads every AF string the same way.
_AF_NUMBER = re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z')


def _rank_af(line, k):
    """How frequent ALT k of this record is, to rank the ALTs of a four-base position:
    the highest of its INFO/AF and every INFO/AF_<x> value (0 when there is none). A
    genotyped panel carries AF; a sites-only merge such as CRISPRme's mega carries only
    per-dataset AF_<x> and AF_max, where AF alone would rank every ALT 0 and leave the
    choice to base order. Ranking only: the dict's AF field is still INFO/AF."""
    vals = []
    for entry in line[7].split(";"):
        key, sep, value = entry.partition("=")
        if sep and (key == "AF" or key.startswith("AF_")):
            v = value.split(",")
            if k < len(v) and _AF_NUMBER.match(v[k]):
                vals.append(float(v[k]))
    return max(vals) if vals else 0.0


def SNPsProcess(line):
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
    # code is that of the site written as one record with every ALT, in any order.
    mask = snp_mask[pos] or _BASE_BIT.get(referenceNucleotide, 0)
    for k, alt in enumerate(altAlleles):
        bit = _BASE_BIT.get(alt, 0) if len(alt) == 1 else 0
        if bit:
            mask |= bit
            alt_af_pos.append(pos)
            alt_af_bit.append(bit)
            alt_af_val.append(_rank_af(line, k))
    snp_mask[pos] = mask
    if mask == 15:
        # all four bases: N would be written, which the search reads as a gap, not as
        # "any base". Resolved once every record is read (resolveFourBasePositions).
        four_base_positions.add(pos)
        return
    letter = _IUPAC_BY_MASK.get(mask)
    if letter is not None:
        genomeList[pos] = letter #inserisco il IUPAC nt nel genoma enriched


def resolveFourBasePositions():
    """Code the positions whose reference base and SNP ALTs cover all four bases.

    The IUPAC letter for A/C/G/T is N, and N is what buildTST and searchTST treat as an
    assembly gap: buildTST discards every target with an N in its guide region and the
    PAM scan matches nothing at an N. Writing it would remove every target through the
    position from the index -- the reference one included. So such a position keeps the
    reference base and its most frequent ALTs, leaving out the rarest ALT by _rank_af
    (a tie leaves out the later of A, C, G, T; an ALT seen on several records counts
    with its highest value). The SNP dict still lists every ALT. On chr22 of HGDP + 1000
    Genomes, 2,448 positions carry three SNP ALTs; the rarest of them has AC >= 41 at 3.
    """
    if not four_base_positions:
        return
    best = {}
    for p, b, v in zip(alt_af_pos, alt_af_bit, alt_af_val):
        if p in four_base_positions:
            key = (p, b)
            if key not in best or v > best[key]:
                best[key] = v
    for pos in sorted(four_base_positions):
        ref_bit = _BASE_BIT.get(genomeStr[pos], 0)
        alts = [b for b in (1, 2, 4, 8) if b != ref_bit]
        drop = min(alts, key=lambda b: (best.get((pos, b), 0.0), -b))
        genomeList[pos] = _IUPAC_BY_MASK[15 & ~drop]
    print("NOTE: %d position(s) of %s carry all four bases between the reference and "
          "the SNP ALTs; N would read as a gap, so each is coded without its rarest ALT "
          "(lowest INFO AF or AF_*, ties to the later of A, C, G, T)"
          % (len(four_base_positions), currentChr))


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


def _af_value(line, pos_AF):
    # Raw AF value string for this record, or "" if unavailable. Handles three
    # cases uniformly: no AF-like field found anywhere in the file (pos_AF is
    # None); this record's INFO field has fewer entries than the position
    # cached from record 1 (variable-length INFO across records); and, as a
    # defensive re-check on every call (not just the initial scan), the cached
    # position pointing at a non-AF key in this particular record. "" is the
    # deliberate placeholder (not "."): it's the only value besides "0" that
    # pandas.to_numeric() parses as NaN instead of raising downstream.
    if pos_AF is None:
        return ""
    splitted = line[7].split(";")
    if pos_AF >= len(splitted):
        return ""
    key, sep, value = splitted[pos_AF].partition("=")
    return value if key == "AF" else ""


def _entry_alt(entry):
    """The ALT of one dict entry "<samples>;<REF>,<ALT>;<rsID>;<AF>", read from the
    right, so a ';' in a sample name cannot move it (the rsID's ';' is already ',')."""
    return entry.rsplit(';', 3)[1].split(',', 1)[1]


def _dict_put(key, entries):
    """Store one record's per-ALT entries under "chr,pos", MERGED with what earlier
    records stored there.

    A position's value is its ALTs' entries joined by '$', the form a record listing
    several ALTs has always produced, so every reader of multi-ALT values reads it. Each
    entry keeps its own record's carriers, genotypes, rsID and AF: the genotypes of a
    split record use 1 for its ALT, as tier0_registry's default alt_index expects.
    Assigning instead let a later record at the position overwrite the earlier ones, so
    a split multi-allelic site kept only its last ALT. An ALT already listed here (the
    same allele on two records) has its entry replaced by the later record's, as the
    whole value used to be; any other ALT is appended, in file order."""
    old = chr_dict_snps.get(key)
    if old is None:
        chr_dict_snps[key] = '$'.join(entries)
        return
    merged = old.split('$')
    alts = [_entry_alt(e) for e in merged]
    for e in entries:
        a = _entry_alt(e)
        if a in alts:
            merged[alts.index(a)] = e
        else:
            merged.append(e)
            alts.append(a)
    chr_dict_snps[key] = '$'.join(merged)


def add_to_dict_snps(line, pos_AF):
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
        af = _af_value(line, pos_AF)
        if len(list_samples) > 0:
            _dict_put(chr_pos_string, [','.join(sorted(list_samples)) + ';' + ','.join(list_chars) + ";" + rsID + ";" + af])
        else:
            _dict_put(chr_pos_string, [';' + ','.join(list_chars) + ";" + rsID + ";" + af]) #None
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
            af = _af_value(line, pos_AF).split(',')

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
            _dict_put(chr_pos_string, final_entry)

def dictSave():
    #os.chdir("./SNPs_genome/")
    with open('my_dict_' + currentChr + '.json', 'w') as f:
        json.dump(chr_dict_snps, f) 


def indel_to_fasta(line, id_indel, pos_AF, start_fake_pos):
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
                af = _af_value(line, pos_AF).split(',')

                start_position = int(line[1])-26
                end_position = int(line[1])+26+len(line[3])
                if start_position < 0:
                    # The window is anchored so the REF base lands at index 25, which a
                    # negative start cannot satisfy: Python would slice from the far end
                    # of the chromosome and hand back an EMPTY window, so the contig and
                    # the coordinates logged for it would disagree. (enricher.cpp throws
                    # instead, so the two ports only agree if both skip.)
                    print(f"WARNING: skipping indel at {currentChr}:{line[1]} -- fewer than 26 bases from the chromosome start")
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
                        allele_fasta = sub_fasta[0:25] + re.sub(line[3], indel, sub_fasta[25:], 1, flags=re.IGNORECASE)

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
    log_indels = []
    id_indel = 1
    start_fake_pos = 0
    if not os.path.isdir(f"fake_{vcfName}_{currentChr}"):
        os.mkdir(f"fake_{vcfName}_{currentChr}")
    output_dir_indels = os.path.abspath(f"fake_{vcfName}_{currentChr}")
    fasta_out = open(f'{output_dir_indels}/fake{currentChr}.fa', 'w+')
    fasta_out.write(f">fake{currentChr}\n")
    list_fasta_indels = []

first_line = True
pos_AF = None  # stays None if no exact "AF=" INFO key is found in the first record
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
    if first_line:
        first_line = False
        splitted = line[7].split(";")
        for pos, ele in enumerate(splitted):
            # exact key match, not a prefix match: a prefix match would also
            # accept e.g. gnomAD-style "AF_afr=", and either garble that value
            # (the "AF="-only extraction below) or -- worse -- lock onto it and
            # skip a later, real "AF=" field in the same record. Since the
            # `break` only fires on a real match, this still keeps scanning
            # past a non-matching population-specific field to find the true
            # "AF=" wherever it appears in the record.
            if ele.split("=", 1)[0] == "AF":
                pos_AF = pos
                break
    if line[6] not in ('PASS', '.'):  # '.' = no filter applied (e.g. HPRC vcfwave VCFs)
        continue
    if sys.argv[4] == 'yes': #if true do all the creation, if false do only genome enrichment with SNPs
        #print(sys.argv[4])
        add_to_dict_snps(line, pos_AF)
        id_indel, start_fake_pos = indel_to_fasta(line, id_indel, pos_AF, start_fake_pos)
    SNPsProcess(line)


resolveFourBasePositions()
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

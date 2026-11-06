# Multi-allelic sites in `add-variants` (CRISPRitz 2.10.0)

What changed in the enricher, what it measured on a real panel, and what CRISPRme must
change, rebuild and test to use it. Written on branch `feat/multiallelic-enrich` (based
on `v2.9.1`) and updated for CRISPRitz 2.10.0, which merged it with
`feat/engine-boundaries` (docs/INDEX_FORMAT.md) and added Q, the four-allele symbol
(section 2.2), the carriers of one allele on several records merged (section 2.3), each
record's AF read by its key (section 2.5) and a reference N kept (section 2.8). The
measurements of section 4 were taken before those four changes, with the four-base rule
that Q replaced. The CRISPRme code cited is the tree at `a5379e0` of the CRISPRme+ 2.8.1
integration branch; line numbers refer to that tree.

**Decided by the owner:** these changes, and the index rebuild they require, go into
CRISPRitz 2.10.0 and the CRISPRme release that pins it, **not** into CRISPRme+ 2.8.1.

## 1. The short answer: were all ALTs of a multi-allelic site enriched?

**No: up to CRISPRitz 2.9.1 only the last one was, wherever a site is split into one record
per ALT — which is how all three published panels are normalised** (measured on chr22 of
each, section 8). They go through `bcftools norm -m -any` (CRISPRme `METHODS.md`,
`docs/HGDP_1KGP_PHASED_INDEX_METHODS.md`), which writes a site with several ALTs as one
record per ALT. The enricher built each position's IUPAC code from the *reference base*
plus the ALTs of the record in hand and overwrote the position, and keyed its SNP
dictionary by `chr,pos` and overwrote that too. So of the records at one position only the
**last in file order** survived, in the enriched genome (hence the search index) and in the
dictionary (hence CRISPRme's Tier-0 registry and Tier-1 genotype store, which are compiled
from it). A record listing several ALTs on one line was coded correctly (except with three
SNP ALTs, section 2.2), but no published panel enriches one: they are normalised first (0 of
1,998,332 chr22 records of HGDP + 1000 Genomes carry a comma in ALT).

**From 2.10.0 every ALT is enriched**, in any record order. A position whose reference base
and ALTs cover all four bases is written `Q`, a fifth symbol the search engine reads as any
base (section 2.2; the IUPAC letter `N` would be read as a gap). The dictionary lists every
ALT, each with its own record's genotypes as they were phased. On HGDP + 1000 Genomes that
includes 672 chr22 sample-positions with two ALTs on one haplotype slot, which CRISPRme's
per-sample analysis turns into a CONFIRMED haplotype holding both (measured through the
CLI): CRISPRme must handle them before it uses a 2.10.0 index (section 6.3).

## 2. What the enricher does now

Both ports, `sourceCode/Python_Scripts/Enrichment/enricher.py` and `enricher.cpp`, byte for
byte (verified, section 5).

### 2.1 The enriched genome

The code at a position is the IUPAC letter of **{reference base} ∪ {every single-base
A/C/G/T ALT of every PASS (or `.`) record at that position}**, whatever the order and
adjacency of those records — the code the same site gets when written as one record listing
every ALT. A record adds its ALTs to what earlier records coded; it no longer starts from
the reference base. The rest of the rule is that of 2.9.1: a member that is not A/C/G/T
(`*`, `N`, a lowercase base) adds nothing, a reference `N` or IUPAC base adds nothing of
its own, and a set of one base leaves the FASTA letter in place.

The running union is kept in a byte per position beside the genome (`snp_mask` /
`snpMask`), because a letter in the genome cannot tell a code the enricher wrote from an
ambiguity code the reference FASTA already holds.

### 2.2 Four bases: Q

The IUPAC letter for {A, C, G, T} is `N`, and `N` is what `buildTST` and `searchTST` treat
as an assembly gap, not as "any base": `buildTST` discards every target with an `N` in its
guide region (`guideRegionHasN`, `mainParallel.cpp`), the PAM scan matches nothing at an `N`
(`pamGenerator.cpp`, `0000`), and the tree has no encoding for it. Writing `N` therefore
removes **every** target through the position from the index, the reference one included.
Demonstrated through the real CLI (`add-variants` → `index-genome` → `search`, NGG, 2
mismatches) on a synthetic chromosome with a perfect target and a reference-`A` site in its
guide carrying ALTs C, G, T:

| enricher | input | code | targets found |
|---|---|---|---|
| 2.9.1 | one record `A>C,G,T` | `N` | **0** (the perfect target is gone) |
| 2.9.1 | three split records | `W` (last ALT only) | 1 |
| `feat/multiallelic-enrich` | either | `V` (A, C, G; T, AF 0.001, left out) | 1, 0 mismatches |
| 2.10.0 | either | `Q` (A, C, G, T) | 1, 0 mismatches |

So such a position is written **`Q`**, a symbol outside IUPAC that `buildTST`, `searchTST`
and `searchBruteForce` read as the mask `1111`: it matches every guide base and every PAM
letter, and is never a mismatch, while `N` stays the gap `0000`. An index file holding `Q` is
flagged in its header (docs/INDEX_FORMAT.md section 4.1, which gives the encoding), so an
engine that cannot read it refuses it. This replaces the four-base rule of the branch,
which kept the reference base and the two most frequent ALTs and left out the rarest (by
`INFO/AF` and every `INFO/AF_<x>`): the owner chose the fifth symbol so that all four
alleles are searchable (CRISPRme plan, decision B1). The left-out ALT's bookkeeping is
gone: nothing is ranked, so no AF is read for the code. The decision is taken as each
record is read (a position that holds four bases holds them whatever comes later), so mode
`no` (genome only) codes exactly what mode `yes` codes. One line on **stdout** counts these
positions per chromosome (`NOTE: <n> position(s) of <chr> carry all four bases ...; each is
coded Q, ...`; stderr is fatal to a CRISPRme stage).

`Q` is also what a record listing three SNP ALTs now gets, where 2.9.1 wrote `N`; a reference
`N` never becomes `Q` (section 2.8); a reference IUPAC base adds no base of its own (2.1),
so `Q` there needs ALTs covering all four bases. A reference FASTA must not hold `Q` (it
would pass a base of its own off as a four-allele position): both ports, and `add-variants`
for a chromosome it copies, refuse one with `ERROR: <chr> holds 'Q' at <pos>: the reference
FASTA must hold A/C/G/T, N and IUPAC codes only (Q is reserved for four-allele SNP
positions)` and exit status 1.

On chr22 of HGDP + 1000 Genomes, 2,448 positions carry three SNP ALTs (section 4): the
four-base rule left out an ALT at each (AC = 1 at 1,193 of them, AC ≥ 41 at 3); `Q` codes
them all. Its cost is the density of windows over such positions (an all-`Q` position
matches any guide base): CRISPRme's ambiguity budget has to count it as four bases.

Measured with 2.10.0 on that chromosome, through `crispritz.py` in the CRISPRme 2.8.1
SIF on ml008 (amd64), the new engine first on `PATH`: `add-variants` (compiled port, dict
mode) on hg38 chr22 and the published `hgdp1kgp.chr22.norm.vcf.gz` (1,998,332 records)
took 979 s and 7.16 GB, wrote nothing to stderr, and coded 2,448 positions `Q`; the
enriched genome holds 11,658,691 `N`, at exactly the reference's positions, and its
sequence is byte-identical to the fifth-symbol prototype's. `index-genome` with NRG and
`-bMax 1` (18 s) stored 18,314 targets holding `Q` and flagged all 3 of its files. Five
guides (CTAACAGTTGCTTTTATCAC, ACTGAAATCTGTAAGCAGGC, CTCTCAGCTGGTACACGGCA,
GAGTCCGAGCAGAAGAAGAA, GGTGAGTGAGTGTGTGCGTG) at 4 mismatches without bulges gave 9,311 rows,
equal to `searchBruteForce`'s on the enriched FASTA, none holding an `N`, every row the
chromosome's bases: 8,954 of them hold `Q`, and 8,951 of those lie in IGL (chr22:22-23 Mb),
where a guide window can hold up to 20 `Q` (276 rows are perfect matches through `Q`).
With one DNA and one RNA bulge: 4,372,248 rows, 4,208,718 through `Q` (326 s). Outside
IGL, `Q` adds 3 rows. These windows are the density the plan's IG guard and ambiguity
budget must bound in CRISPRme (a `Q` is two bits of ambiguity, as many as two IUPAC
letters of two bases).

### 2.3 The SNP dictionary (`my_dict_<chr>.json`)

The value of `"<chr>,<pos>"` lists **every ALT** of the position: one entry per ALT joined
by `$`, the form a record listing several ALTs has always produced, so every reader of
multi-ALT values reads it unchanged. Each entry is `<sample:GT,...>;<REF>,<ALT>;<rsID>;<AF>`
built from **its own record**, exactly as 2.9.1 built a lone record's value: its carriers
(sorted), their genotypes verbatim (for a split record `1` is that entry's ALT, which is
`tier0_registry`'s default `alt_index`), its rsID (`;` normalised to `,`) and its own AF.
Entries are in file order. A second record of an allele already listed (the same `REF>ALT`
twice: 362 genome-wide in the published panels, 351 pure duplicates) **merges into that
allele's entry**, in place: the entry lists the carriers of both records; a sample listed
by both, phased alike (as many alleles, every separator `|`), carries the ALT on every slot
where either record puts it (`0|1` and `1|0` give `1|1`), and otherwise, as unphased slots
cannot be matched up, keeps the genotype with more copies of it (on a tie, the earlier
record's); and the entry takes the rsID and AF of the more frequent record (the higher
`INFO/AF`, a missing or non-numeric one counting 0; on a tie, the earlier). A merged entry's
genotypes are in split form, `1` being its ALT and `0` the reference or another ALT, also
where one of the records listed several ALTs (its `1|2`, for its second ALT, reads `0|1`
there), so they never mix two records' allele numbering. Up to the branch the later
record's entry replaced the earlier one and its carriers were lost (owner decision B2:
merge; section 7). CRISPRme's Tier-1 compiler refuses a repeated `(pos, alt)`
(`tier1_genotypes.py:283-287`), so the dictionary never holds one. Keys and key order are
unchanged. A key holding one record is byte-identical to 2.9.1, but for its AF (section
2.5).

The genotypes are copied verbatim, phase included, so split records that were phased
independently come through as they were phased, even where that puts two ALTs of one
position on one haplotype slot (`0|1` under G and `0|1` under T): 672 sample-positions on
chr22 of HGDP + 1000 Genomes (section 4). The enricher neither resolves nor drops such a
conflict (unit test `test_dict_passes_conflicting_phase_through`, `verify_enricher.sh` case
30); whatever builds haplotypes from the dictionary has to. For CRISPRme see section 6.3.
Merging the carriers of an allele on two records can make such conflicts too (on HPRC
chr22, 37 sample-positions at 38,152,965; section 7, decision 2): CRISPRme resolves them
with the same rule (owner decision: two ALTs on one haplotype slot are PUTATIVE).

### 2.4 The indel log and the fake-indel genome

Unchanged by the merge and by `Q`: `indel_to_fasta` cuts every contig from the plain
reference (`genomeStr`), never from the enriched genome, so a merged SNP code next to or at
an indel does not reach a contig, its log row or its FAKEPOS (unit test
`test_indel_outputs_unmoved_by_split_snps`; `verify_enricher.sh` cases 19, 21, 28;
byte-identical on chr22 of all three panels, section 4, and with 2.10.0 too). The log's AF
column changes where a record lists `AF` at another place in its INFO than the file's first
record, as each record's AF is now read by its key (section 2.5): HPRC's chrY log gains the
AF of 91,950 of its 102,665 rows, every other column unchanged, so it is rebuilt (section
8). The multi-allelic indel fix of 2.9.1 (each allele's contig from the reference window) is
untouched and still covered by case 7. CRISPRme's overlay (`overlay_indel_snps.py`) copies
the enriched genome's codes onto the contigs' flanks afterwards, but only those its
`_IUPAC_AMBIG` lists, which up to CRISPRme 2.8.1 has no `Q`: it must learn `Q`, or every
flank base at a four-allele position keeps its reference base (section 6.4).

### 2.5 FILTER and AF act per ALT, and AF is read from each record

A record whose FILTER is neither `PASS` nor `.` contributes nothing, to the code or the
dictionary, and the PASS records at its position are still merged (test
`test_non_pass_record_contributes_nothing`, case 23). Each ALT keeps its own AF in the
dictionary. Each record's AF is its `INFO/AF` value, found by its key wherever it stands in
that record's INFO (an exact `AF`: an `AF_afr` or `AF_max` is not it); a record without one
gets an empty AF, in the dictionary and the indel log. Up to 2.9.1, in both ports, the AF
field's place was taken from the file's first record, so every record whose INFO lists its
fields in another order or number got none: the retired 1000G2021_HGDP index (a merged VCF)
lost the AF of 224,981 of 232,384 chr22 indels that way, and HPRC chrY 91,950 of its
102,665 indel AFs (unit test `test_af_read_per_record_by_key`, `verify_enricher.sh` case
31). A minor-allele-frequency filter applied before enrichment acts per ALT only on split
records — measured with bcftools 1.24 on a site `G>A,C,T` with AF 0.3, 0.0005, 0.02 and
`bcftools view -i 'INFO/AF>0.001'`:

| input to `add-variants` | 2.9.1 code / dict ALTs | 2.10.0 code / dict ALTs |
|---|---|---|
| split, unfiltered | `K` / T | `Q` / A, C, T |
| split, then filtered (C dropped) | `K` / T | `D` / A, T |
| one record, filtered (kept whole: one ALT passes) | **`N`** / A, C, T | `Q` / A, C, T |

So filter **after** splitting, as the published panels are built.

### 2.6 Inputs with one SNP record per position

What sections 2.1 and 2.3 change needs two SNP records at a position, so it changes nothing
where no position holds more than one: on the 17-case `verify_enricher.sh` fixture of 2.9.1
with its two split sites reduced to one record each, both ports of 2.9.1 and of the branch
gave identical md5 for all four outputs and the run log. Such an input is still not enriched
by 2.10.0 as by 2.9.1 wherever

* a record lists three SNP ALTs: `Q`, where 2.9.1 wrote `N` (2.2);
* a record lists two or three ALTs at a reference `N`: the `N` stays, where 2.9.1 coded it
  (2.8);
* a record lists `AF` at another place in its INFO than the file's first record: it gets
  its AF, in the dictionary and the indel log (2.5);
* a SNP lies outside the chromosome: an error (2.7), as is a reference holding `Q` (2.2).

Measured with 2.10.0 against 2.9.1 (the CRISPRme 2.8.1 SIF's) on ml008. On the reduced
fixture both ports of each give the same dictionary, indel log, fake-indel genome and run
log, and genomes that differ at one base, 1,618, where case 17 lists `N>A,C` at a reference
`N` (`M` in 2.9.1, `N` now). Through `add-variants` (compiled ports) on hg38 chr22 with the
raw HPRC vcfwave VCF (`ALL.chr22.HPRC.vcfwave.vcf.gz` of `vcfs/hg38_HPRC/`, 623,085
records, which keep a site's ALTs on one record: 529,585 positions with a SNP record, none
with two; no record listing two or three ALTs at a reference `N`; `AF` at the same place in
every record's INFO), the dictionary, indel log and fake-indel FASTA are byte-identical, and
the genome differs at 61 bases, those of the 61 records listing three SNP ALTs: `N` in 2.9.1
(11,658,752 `N` against the reference's 11,658,691), `Q` now. So an enrichment of a VCF that
keeps a site's ALTs on one record changes too wherever one of the above occurs, as on the
raw HPRC vcfwave VCFs: rebuild it.

### 2.7 A SNP outside the chromosome

A SNP record whose POS is 0 or past the last base (a VCF of another assembly) now stops the
enricher with `ERROR: <chr>:<pos> is outside the chromosome (<n> bp)` on stderr and exit
status 1, in both ports; `add-variants` then fails naming the chromosome. In 2.9.1 the
Python port raised an `IndexError` past the end and coded the **last** base for POS 0, and
the C++ port read and wrote out of bounds; the new per-position mask would have made that
write unconditional.

### 2.8 A reference N stays N

An `N` in the reference is a gap, and no VCF record makes it a base: the enriched genome
holds an `N` exactly where the reference does, so its N count is the reference's (unit test
`test_enriched_genome_keeps_the_references_N`, the reconstruction check of
`verify_enricher.sh`, and a rebuild gate). Up to 2.9.1 a record listing two or three ALTs at
a reference `N` (`N>A,C`) was coded over it (`M`), making a gap searchable, and a record
listing three SNP ALTs wrote `N` over a base. Split records at an `N` were not coded over it,
as 2.9.1 coded each from the reference base alone, but merging them (section 2.1) would have
been, as on the branch (`f266eac`). The dictionary still lists a record at an `N`.

## 3. Cost

One byte per base for the union, in both ports; the dictionary merge runs only at a key
already present. Measured on chr22 with the four-base rule, which kept ~17 bytes per
single-base ALT more for its ranking (section 4); `Q` keeps none.

## 4. Measured on real chromosomes: chr22 of the three panels

Inputs: the published panel's raw VCF (`hgdp1kgp.chr22.norm.vcf.gz`, 1,998,332 records,
4,091 samples, FILTER `.` throughout, 0 commas in ALT, 0 `*` ALTs) and the hg38 chr22 FASTA
of the published install (REF equals the FASTA base at all 1,757,775 SNP-bearing records).
Each run is the real `crispritz.py add-variants <vcfdir> <genomedir> true`, on ml008 in the
CRISPRme+ 2.8.1 SIF: **before** = the shipped CRISPRitz 2.9.1 (its compiled `enricher`),
**after** = this branch, once with the compiled port and once with `enricher.py`. The
expected outputs come from an independent path: a bcftools 1.24 query of the same VCF
(site fields of every record; genotypes of every record at the multi-SNP positions),
compared by a separate script.

Census (bcftools site query): 1,711,938 positions carry a SNP-bearing record; **43,389**
carry two or more (40,941 two, 2,448 three); **45,837** SNP records (2.6 %) are not the last
at their position — 5,104 of them with AC ≥ 41 and 1,566 with AC ≥ 818 (10 % of 8,182
haplotypes). No position repeats an ALT. At 288 positions the SNP records are not
adjacent; 1,289 positions also carry an indel record. The published `registry_hg38_HGDP1kGP`
(the index-only install's, compiled from a pre-2.10.0 dictionary) holds 1,711,938 positions on
chr22 — exactly the positions with a SNP record — none with two ALTs.

Phase of the split records (counted in the branch's dictionary, and again from a bcftools query
of the VCF's genotypes; both give the same numbers). At the 43,389 positions, 55,826
sample-positions carry two or more different ALTs: 55,154 on different haplotype slots, as a
`1/2` heterozygote should be, and **672 with two ALTs on one slot**, at **584 positions**, in
405 samples. The panel was phased (SHAPEIT5) as biallelic records — its
`bcftools norm -m -any` step found nothing left to split (CRISPRme
`docs/HGDP_1KGP_PHASED_INDEX_METHODS.md`, section 7.1) — each record a variant of its own,
and nothing in that keeps a haplotype to one ALT per position: e.g. 10,661,629 NA18869,
`0|1` under G and `0|1` under T; 10,683,738 HG02804, `1|0` under C and `1|0` under T. 91 of
the 672 carry more ALT copies than the two slots hold (`1|1` and `0|1`); 220 sit at a
position with three SNP ALTs. Counting each extra ALT on a slot, there are 746 slot
conflicts. None of the 55,826 has an unphased genotype. At 231 (position, ALT) pairs every
carrier of the ALT is one of these samples (e.g. 10,723,367 G, carried only by HG01891). Up
to 2.9.1 the dictionary held one ALT per position, so nothing downstream could meet a
conflict; section 6.3 says what CRISPRme must do with them.

Before vs after (compiled ports; the Python port, below, is byte-identical to the compiled
one):

| output | result |
|---|---|
| enriched genome | **43,389 bases change, exactly the multi-SNP positions**; every base equals the rule: 2.9.1 = reference + last record's ALT, the branch = the merged rule of 2.1 with the four-base rule, at all 1,770,802 positions with a REF-length-1 record (0 mismatches each); no base outside SNP positions differs from the reference in either; `N` written for a variant: 0 and 0 |
| codes | all 2-base → 3-base: Y→H 9,349, S→V 8,623, K→D 8,262, Y→B 6,943, K→B 3,838, W→D 2,175, R→V 2,079, W→H 876, K→V 526, W→V 317, K→H 210, Y→V 191 |
| non-last ALTs | all 45,837 absent from the 2.9.1 code; 44,633 present in the new code; 1,204 are the rarest ALT of a four-base position |
| four-base positions | 2,448 (the stdout note says so; 2.10.0 codes each `Q`). At 1,244 the rarest ALT is the last record, which 2.9.1 coded and the branch left out for the two others (AC = 1: 626, 2-9: 571, 10-40: 46, ≥ 41: 1); at 1,204 it is an earlier record (AC = 1: 567, 2-9: 574, 10-40: 61, ≥ 41: 2) |
| SNP dict | same 1,711,938 keys in the same order; **43,389 values change, exactly the multi-SNP positions**; each equals the bcftools reconstruction — after: every record's entry in file order (43,389/43,389), before: the last record's entry (43,389/43,389); 45,837 entries added |
| indel log, fake-indel FASTA, `change_version.txt` | byte-identical (md5 `4197a2d0…`, `34ed6062…`, `a824c4b2…`) |
| CRISPRme overlay (`overlay_indel_snps.overlay_fake_chromosome`, a5379e0) | same length (13,693,118), 686,631 flank bases coded in each; **50,253 differ**, every one at a multi-SNP position and equal to its new code |
| cost, compiled port | 986 s → 1,004 s wall; peak RSS 6.72 → 6.86 GB |

**The two ports at full scale:** `enricher.py` through the same CLI (1,789 s, peak RSS
3.66 GB) gives the same five files, md5 for md5 (enriched genome `21903a97…`, dict
`16fd0e8c…`, log, fake-indel FASTA, `change_version.txt`), and the same run log. The branch's
final code (commit `8e36f57`, which added the `AF_*` ranking) gives the same five md5s again
(this panel's INFO holds only `AC`, `AN`, `AF`), as do the CLI demonstrations of 2.2 and 2.5.

**The other two panels**, chr22, the same way (2.9.1 vs the final branch, through
`add-variants`, checked by a second script that rebuilds the expected code and dict values —
split and multi-ALT records, the AF/AF_* ranking — from a bcftools query of the VCF):

| panel (chr22 VCF) | multi-record positions | four-base | bases changed (= the expected set) | dict values changed | multi-record dict values equal to the reconstruction (after, before) | log, fake-indel FASTA |
|---|---|---|---|---|---|---|
| HPRC, `ALL.chr22.HPRC.vcfwave.vcf.gz` normalised as `hprc_build.sh` does | 4,765 | 92 | 4,761 | 4,761 | 4,765, 4,765 | byte-identical |
| mega, `mega.chr22.afmax.vcf.gz` (sites-only) | 16,103 | 1,370 | 16,103 | 16,103 | 16,103, 16,103 | byte-identical |

In both, every base of both genomes matches its rule (0 mismatches) and no variant is coded
`N`. HPRC's four unchanged positions hold one ALT on two records (section 7, decision 2).
HPRC's split records are phased consistently: 7,666 sample-positions carry two ALTs at one of
its 4,761 positions with two or more SNP ALTs in the dictionary, none on one slot (counted as
for HGDP + 1000 Genomes above, both ways), with the later record of a repeated allele kept as
on the branch; 2.10.0 merges the carriers of such records (section 2.3), which puts two ALTs
on one slot for 37 sample-positions, all at 38,152,965. Mega is sites-only: no genotypes.

## 5. Tests on the branch

* `test/test_enricher_multiallelic.py` — 24 unit tests, each run against **both** ports
  (48 in all), with exact expected outputs: two and three ALTs in every order, split ==
  joined, four bases coded `Q` (every order, whatever the AF fields, a record listing three
  ALTs, an ALT on two records), a reference holding `Q` refused, reference `N` and IUPAC
  bases, the N-count invariant (an `N` exactly where the reference has one, in both modes),
  mode `no` == mode `yes`, the dictionary entries (carriers, rsIDs, AFs, haploid / unphased /
  missing / `GT:DP` genotypes, records without carriers, two or three ALTs on one haplotype
  slot passed through verbatim), a non-PASS record, the same allele on two split records
  and on a split record and a record listing several (merged carriers), each record's AF by
  its key, a split record beside a multi-ALT record, `*`, non-adjacent records, the indel
  outputs, a position with one record, and a SNP outside the chromosome. Against the
  branch's last commit (`f266eac`, the four-base rule) 26 of the 48 fail, against 2.9.1 44
  (the 4 that pass, mode `no` == mode `yes` and a position with one record, pin what does
  not change).
* `verify_enricher.sh` — grows from 17 to 33 cases (cases 18-30: split sites in every
  order, the same sites as one record, four bases, non-adjacent records, genotype shapes,
  SNP + indel at one position and in an indel's window, a split record beside a multi-ALT
  record, a non-PASS record, the same allele twice, AF spellings, `AF_<x>` fields, a
  reference IUPAC base and `N`, `*`, the chromosome edges, two and three ALTs on one
  haplotype slot; 31-33: AF read by its key, the same allele on two records merged, ALTs on
  a reference `N`); both ports byte-identical, with the file set and the run log, and
  checked against what they must be: the enriched genome against an independent
  reconstruction from the FASTA and the VCF (with the N-count invariant), the dictionary
  and the indel log at cases 24 and 31-33, and a reference holding `q` refused by both
  ports with one message. Earlier, on the branch, also checked against an independent
  oracle (91 coded positions, 81 dictionary keys, 0 mismatches) and an independent carrier
  check at the 62 positions whose SNP records each list one ALT, none repeated (each entry
  lists exactly its record's carriers, genotypes verbatim; 0 mismatches).
* `test/run_four_allele_test.sh` — `Q` end to end through `crispritz.py` (`add-variants`
  with both ports, `index-genome`, `search`), against an independent oracle and
  `searchBruteForce`; docs/INDEX_FORMAT.md section 11.
* `test/run_add_variants_status_test.sh` — with the shipped enricher, whatever is on
  `PATH`: a chromosome `add-variants` copies (no VCF) holding `q` is refused before any
  enrichment, a README, `.DS_Store`, `.fai` and macOS `._<name>.fa` beside the chromosomes
  are ignored (exit 0, nothing on stderr), and a folder holding no FASTA is refused.
* The existing scripts pass unchanged on Linux (g++ 14): `run_minimal_test.sh` (NGG, WTN,
  WTN under AddressSanitizer), `run_golden_test.sh`, `run_pam_filter_test.sh` (31
  assertions), `run_flank_n_test.sh`, `run_short_contig_test.sh`.

## 6. What CRISPRme must change

Nothing in CRISPRme's file formats has to change: the Tier-0 registry and the Tier-1 store
already hold one record per `(pos, alt)` and several per position, and every dictionary
reader already splits a value on `$`. What changes is **content** (positions that held one
ALT now hold two or three, and the enriched genome holds `Q`), so the data must be rebuilt
(section 8), a few places that assumed one ALT per position, or that were only ever
exercised with one, must be fixed, and every table that decodes IUPAC letters must read `Q`
as any base (docs/INDEX_FORMAT.md section 8, item 8; in the SNP path, 6.3, and the overlay,
6.4).

### 6.1 Pin, provenance, and a guard against mixing

* Pin CRISPRitz `2.10.0`: `Dockerfile:31` (`ARG crispritz_ref=v2.9.1`), `environment.yml`
  and `install_from_source.sh` (which name 2.9.1), the Bioconda `crisprme` recipe.
* The enriched genome, its TST index, the fake-indel genome and its `_INDELS` index, the
  per-sample dictionaries and the Tier-0/Tier-1 stores must all come from **one**
  enrichment. Mixed, the code and the genotype source disagree at ~2.5 % of SNP positions:
  an old dictionary under a new genome has no carriers for a coded ALT (rows needing it are
  not emitted), a new dictionary under an old genome offers ALTs the search never matched. Stamp the enrichment in the tier manifests
  (`registry_<vcf>/reg_<chrom>.idx`, whose manifest `tier0_registry._write_registry`
  already stamps with `data_type`/`phased`, line 1158, and `genotypes_<vcf>/gt_<chrom>.idx`,
  written by `tier1_genotypes.compile_genotypes`) and in the index's `manifest.json`, e.g.
  `"enricher": {"crispritz": "2.10.0", "multiallelic": "merged"}`, and refuse (or warn on)
  a search whose index and tiers disagree. `add-variants` is called at `crisprme.py:4150`;
  the tiers are emitted at `crisprme.py:4326-4475`
  (`build_dictless_tiers.emit_dictless_tiers_guarded`).

### 6.2 Dict-less tiers and registry

* `PostProcess/build_dictless_tiers.py` (`emit_dictless_tiers`, lines 165-252) →
  `tier0_compile.compile_from_dict` (lines 268-350) and
  `tier1_genotypes.compile_genotypes_from_dict`: **no code change**. Each `$` entry becomes
  one `(pos, alt)` record (`tier0_compile._emit_key`, lines 247-262); counts use
  `alt_index="1"` (`tier0_compile.py:269`, `tier0_registry._classify`, lines 157-174), which is
  right for split-record entries. Verified by compiling a merged value with CRISPRme's own
  compilers: two registry records at the position, AC/AN per ALT exact, Tier-1 carriers per
  ALT with their genotypes. **Rebuild** them from the new dictionaries.
* Pre-existing, unrelated to split records: an entry from a record listing several ALTs
  keeps that record's genotypes (`1|2`), which `alt_index="1"` counts as 0 copies of the
  second ALT. Not reachable on a published panel; worth a test if such VCFs are accepted.
* The mega registry is compiled from INFO/AF straight from the VCF
  (`compile_registry_from_info_af`), so it already lists every ALT — 20,517 chr22 positions
  with two or more SNP ALTs in the published `registry_hg38_mega` — while the mega enriched
  genome coded one of them. It needs no rebuild; the mega genome and indexes do.
* Side finding, pre-existing and separate: that registry takes **every** record, whatever its
  FILTER (815,676 chr22 SNP positions, 20,517 of them multi-ALT: exactly the counts over all
  records of `mega.chr22.afmax.vcf.gz`), while the enricher codes only `PASS`/`.` records
  (742,436 and 16,103). So the registry lists ALTs of non-PASS records that the genome never
  codes; CRISPRme's registry-only path offers them as candidates. Apply the enricher's FILTER
  gate in `build_mega_registry.py`, or decide that non-PASS mega records should be coded.

### 6.3 The SNP path (`new_simple_analysis.py`)

* `retrieveFromDict` (lines 279-305) → `simple_analysis_registry.retrieve_5tuple` (lines
  103-226), `legacy_5tuple_from_entry` (25-54), `dict_alt_to_samples` (76-100): already
  return one element per ALT. No change.
* Legacy dictionary path, `iupac_decomposition` (lines 898-): each ALT is its own seed
  (`totalDict[...][(pos_c, elem)]`), and the lattice only combines strictly increasing
  columns (`if newkey[-2] > key[-2]`, line 1198), so two ALTs of one position never share a
  haplotype. The per-haplotype test `allele != "0"` (lines 947-952) is right for split-record
  entries (a slot holding the other ALT reads `0`). No change for split records (but see `Q`,
  below); note it is **wrong for entries of a record listing several ALTs** (`1|2` puts the
  sample on both slots of both ALTs) — pre-existing, not reachable on the published panels;
  the observed path below resolves it per record and could be reused here. This path (and
  its population-level greedy, `by_pos`, from line 1027) substitutes **every** ALT the
  dictionary lists, without checking the target's IUPAC column (with the four-base rule it
  also tried the ALT the code left out; with `Q` the column holds every ALT). Filtering by
  the column, as the observed path does, would make the two paths agree.
  A phased sample whose slot holds two ALTs of one column (section 4) is put on that slot's
  seed of both ALTs (lines 944-962), and the peel (lines 1242-1267) leaves it on both
  maximal combinations, so it is listed as a carrier of both ALTs' haplotypes, each row
  naming one ALT per column. This path labels no row CONFIRMED (the only phase it gives is
  PUTATIVE, to the SNP+SNP companion, line 1376), so listing the sample under both is the
  conservative reading: no change.
* **Must change — `Q`, in both paths.** A target goes to `iupac_decomposition` only if it
  holds a letter of `iupac_nucleotides` (line 2068; the set at line 149), and a column is
  decomposed only if its letter is in `iupac_code` (line 927 in the legacy path, line 638 in
  the observed path below; the table at lines 183-205). Neither has `Q`, nor has
  `iupac_code_set` (lines 150-180; the observed path's column filter, line 655), so a target
  whose only variant column is `Q` is taken for a reference target, and at a `Q` column of
  any other target no ALT is substituted: none of a four-allele position's three ALTs
  reaches a row. Add `Q` and `q` to the three as any base, `{A, C, G, T}`
  (docs/INDEX_FORMAT.md section 8, item 8).
* Observed-haplotype path, `_iupac_decomposition_observed` (lines 607-), taken whenever a
  Tier-1 store is present (line 904); it enumerates haplotypes under the per-sample analysis,
  the default on a phased genotyped index such as HGDP + 1000 Genomes (`crisprme.py`, lines
  2134-2148), and returns early with greedy representatives under the population-level one
  (line 699). The alt index is resolved per record from its carriers
  (line 675, `observed_haplotypes._alt_index_for_record`) — `1` for split-record entries —
  and the target's own IUPAC column filters the ALTs (#139, line 655): at a four-allele
  position that column is `Q`, which keeps every ALT once the tables read it as all four
  bases (the `Q` item above); the four-base rule's code dropped the left-out ALT there. Its
  CONFIRMED sets are right **only where a sample's split records put its ALTs on different
  slots** (measured with CRISPRme's own module: `1|0` under A and `0|1` under C gives the
  two sets {A} and {C}). Records phased independently do not always do so, and the default
  panel has 672 such sample-positions on chr22 (section 4):
* **Must change — `observed_haplotypes._per_sample_variant_sets`, CONFIRMED branch (lines
  200-218).** It groups a phased sample's ALTs by slot (`slots.setdefault(h,
  set()).add((pos_c, alt))`, line 208), so a slot holding two ALTs of one column becomes
  **one CONFIRMED set with both**: a row naming two alleles at one position, labelled
  confirmed cis, while the sample is a carrier of neither single-ALT set. Not reachable with
  a 2.9.1 dictionary (one ALT per position); reachable with 2.10.0 on HGDP + 1000 Genomes at
  584 of the 43,389 chr22 positions with two or more SNP ALTs (1.3 %), in 405 samples, and on
  HPRC chr22 at one position (38,152,965, 37 samples), where the merged carriers of one
  allele's two records put two ALTs on one slot (section 2.3; the records' own phase puts
  none, section 4). Measured on the real data:
  - with CRISPRme's own code, as a measurement harness (not the CLI): the 30 dictionary
    keys of chr22:10,661,500-10,661,800 from this branch's chr22 enrichment, compiled by
    `build_dictless_tiers.emit_dictless_tiers` into a Tier-0 registry and Tier-1 store, read
    back the way `_iupac_decomposition_observed` reads them (lines 636-685) for the target of
    guide `TTATTCCAGACAGATTATGC` + `TAG` at chr22:10,661,610 (column 19 = 10,661,629, coded
    `B`; G has 2,133 carriers, T has 2, NA18868 `1|0` and NA18869 `0|1`), then
    `enumerate_observed_haplotypes` (line 828): five sets, one of them `{(19,G),(19,T)}`
    CONFIRMED with NA18869 its only carrier, its sequence holding T, its info both variants
    (`chr22_10661629_C_G,chr22_10661629_C_T`); NA18869 is in neither {G} nor {T};
  - through the real CLI, as a user would meet it: `crisprme.py build-index-only` (this
    branch's compiled enricher in place of the image's, so `add-variants` is the branch's) on hg38
    chr22:1-10,700,000 with the 4,433 HGDP + 1000 Genomes records it holds (every sample),
    the index-only layout, then `complete-search --per-sample` with that guide (NRG, 4
    mismatches, 1/1 bulges); CRISPRme a5379e0 in the rc281 image; exit 0, stderr empty. The
    per-haplotype report (`observed_haplotypes_raw.tsv`) and
    `all_results_with_alternative_alignments.tsv` hold a row at 10,661,609 (+, no bulge, 1
    mismatch) whose variant columns read `chr22_10661629_C_G,chr22_10661629_C_T` and
    `20C>G,20C>T` (MAF `0.000244439,0.33757`, the reverse of the variants' order: G's AF is
    0.33757), samples `NA18869`; the G row lists 2,126
    samples, not NA18869. `phase_confirmation.tsv` labels that haplotype CONFIRMED in all 190
    alignments that carry it (with and without bulges), and `snp_snp_cooc.tsv` reports the
    same 190 as SNP+SNP co-occurrences, CONFIRMED, carrier NA18869: two SNPs "in cis" at one
    position. (`integrated_results.tsv` shows this cluster as its reference perfect match,
    so none of its variant rows reaches the main table in this run.) `--population-level`
    on the same index (one candidate per column, `twopass_emit`): exit 0, no row names two
    ALTs.

  **Fix:** never build a set holding two ALTs of one column. Recommended: a sample whose
  slot holds two ALTs of one column is not phased in that window, so send it to the PUTATIVE
  branch (below), as the trigger at lines 193-200 already does for a sample with an
  unphased genotype or with columns in two phase sets — a third condition there. That is the
  rule the indel path of `r281f-indel` already applies (6.4). The alternative, one set per
  ALT (the slot's other variants with G, and with T), keeps every row well-formed but must
  not label them CONFIRMED either: at least one of the two calls is on the wrong slot, and in
  a `0|1` + `0|1` conflict the other ALT belongs on the other slot, a combination one set per
  ALT does not offer. Either way the sample stays a carrier of {G} and of {T}. Consequence of
  any handling that does not call such a sample CONFIRMED: the enumerator labels a set
  CONFIRMED only if every carrier reached it CONFIRMED (lines 283-303), so a set the sample
  shares turns PUTATIVE for all its carriers — at 10,661,629, {G} with 2,126 carriers. Keeping
  the others CONFIRMED needs a phase state per carrier rather than per set (section 7,
  decision 4).
* **Must change — the PUTATIVE branch of the same function (lines 223-249).** An unphased
  carrier of two ALTs of one position (a `1/2` heterozygote, split into `1/0` under one ALT
  and `0/1` under the other) gets every subset of its carried variants as a candidate,
  **including the set holding both ALTs of that column**: measured, the enumerator emits
  `{(3,A),(3,C)}` with a sequence holding C and info listing both variants. Unreachable while
  the dictionary held one ALT per position. Reachable with 2.10.0 on any panel with unphased
  genotypes, and, once the CONFIRMED fix above is in, on phased panels too: each sample it
  downgrades takes this branch (on HGDP + 1000 Genomes chr22, none of the 55,826 multi-ALT
  sample-positions is unphased, so there it is the 672 conflicts that reach it). Fix: build
  the combinations with at most one ALT per column (`itertools.product` over each carried
  column's {absent, ALT₁, ALT₂, …} minus the all-absent choice), or drop combinations whose
  columns repeat. The fallback for a sample carrying more than `max_putative` variants
  (`IUPAC_CAP`, default 10; line 243) emits the union of everything carried, both ALTs of a
  column included, and needs the same rule (e.g. one union per choice of ALT at such a
  column).
* **Should change — `population_summary`** (the companion's frequencies). A single-variant
  row (k = 1) is read from the registry, exact per ALT. A combination row (k ≥ 2) on a phased
  panel counts a carrier's slot when it holds the ALT at every row position
  (`_phased_decide`, lines 648-676); `summarize` (lines 462-635) reads only the row's
  `(pos, alt)` pairs, so a slot that also holds another ALT of a row position is counted, in
  the row of each ALT. Once the rows follow the CONFIRMED fix, decide such a carrier with
  `_unphased_decide` (a carrier, no haplotype count), finding it through the position's other
  ALTs (`tier0_reader.alts_at`) and their Tier-1 carriers. `twopass_emit` (the
  population-level mode: one candidate per column, rows PUTATIVE): no change.
* Performance: `countIUPAC` and `CRISPRME_IUPAC_CAP` count **positions**, but each position
  now contributes one seed per ALT to the legacy 2^k lattice. At ~2.5 % of SNP positions
  the width grows from 1 to 2-3 seeds; re-measure a dense-window run before release.

### 6.4 The indel path

* **At `a5379e0`, must change** (`analisi_indels_NNN.py`, SNP+indel co-occurrence):
  `_snp_at` (lines 965-981) takes `alt = alts[0]  # biallelic SNP position` — the
  alphabetically first ALT of the registry — so where a target uses another ALT at an
  overlaid column, that SNP is not counted as used: the indel row is still emitted, but the
  SNP+indel co-occurrence (companion TSV, CONFIRMED-cis variant columns and carriers) misses
  it; the AN of the joint AF (line 995,
  `_reg.alts_at(...)[0]`) and the PUTATIVE marginal AF (line 1054) read the first ALT too.
  `indel_snp_cis.used_snps_for_target` (lines 223-259) accepts one `(ref, alt)` per column
  and `forward_snp_tokens` (line 144) picks `alts[0]` again.
* **The `r281f-indel` branch already does this**: `_snp_record` (its
  `analisi_indels_NNN.py:876-895`) returns `{alt: {sample: GT}}` for every ALT of the
  position (Tier-1 or the dictionary), the names and AFs are looked up per ALT
  (`_snp_name(pos1, alt)`, `_snp_global(pos1, alt)`), `SnpColumn.alts` holds every ALT, and
  `indel_snp_cis.biallelic_site` (commit `cd5bb80`) reads a record listing several ALTs
  against each ALT; a split position "comes back untouched". Land it before or with the
  2.10.0 pin. Left over there: `forward_snp_tokens` still picks `alts[0]` but has no caller —
  delete it. The column test `if col.ref not in col.code` stays true once `Q` reads as all
  four bases: every code the enricher writes holds the reference base. Its `resolve_carriers` (`indel_snp_cis.py`, lines 534-665
  at `r281f-indel` e61a761) already refuses a slot holding two ALTs of one column as a
  haplotype (line 629: `slots = None  # two alts on one slot: not a haplotype`) and resolves
  that sample PUTATIVE on the alleles it carries (lines 643-665), which is what 6.3 asks of
  the SNP path.
* **Must change — `overlay_indel_snps.py`** (the same in CRISPRme 2.8.1, `0b8943e`). It
  copies a flank base's code from the enriched genome only if the code is in
  `_IUPAC_AMBIG` (line 34: `RYSWKMBDHV`, upper and lower case), which holds `B`/`D`/`H`/`V`
  (the enricher never writes `N` for a variant) but not `Q`. So a flank base at a
  four-allele position keeps its reference base, and the indel search finds a target that
  needs an indel and one of that position's three ALTs only with a mismatch there, if at
  all. Add `Q` and `q` to `_IUPAC_AMBIG` (docs/INDEX_FORMAT.md section 8, item 8). Measured
  with the script itself, as shipped and with `Q` added, on the 2.10.0 chr22 enrichments
  (in the CRISPRme 2.8.1 SIF on ml008): on HGDP + 1000 Genomes (2,448 `Q`) it codes 671,160
  flank bases, 15,471 fewer than with `Q` (686,631, as many as it codes on the 2.9.1
  enrichment, section 4), and those 15,471 keep their reference base where they should read
  `Q`, in 2,027 of the 240,557 contigs; on HPRC, normalised as `hprc_build.sh` does (92
  `Q`), 613,366 against 624,345: 10,979 bases in 3,344 of 195,813 contigs. The overlaid
  fake-indel genome and its `_INDELS` index must be **rebuilt** (section 4: measured which
  flank bases change).
* `build_indel_genotypes.py` (keyed `<pos>_<REF>_<ALT>`) and the log-driven carriers: not
  affected.

### 6.5 The report

No format change: a row's variant columns list the ALT each haplotype holds, and two ALTs of
one position appear in two rows, as they would for a record listing several ALTs
(`resultIntegrator.py` splits the tokens per row, lines 455-520) — once 6.3 is fixed. Until
then a 2.10.0 index puts rows naming two ALTs of one position into the per-sample report and
its companions (measured, 6.3). The SNP+SNP companion took such a row for a co-occurrence
because `_record_snp_snp_cooc` (`new_simple_analysis.py`, lines 387-391) counts distinct
variant tokens (`chr22_10661629_C_G`, `chr22_10661629_C_T`), not distinct positions; it
should compare positions. Expect more variant rows and more carriers at the affected
positions; the changelog and `METHODS.md` must say that earlier releases reported only the
last ALT of a split site, and what `Q` is (a target may show it in its DNA column).

### 6.6 Tests CRISPRme should add

* A Tier-0/Tier-1 compile of a merged value (two ALTs, carriers on opposite slots): per-ALT
  AC/AN and carriers (the measurement in 6.2, as a unit test).
* The observed CONFIRMED fix (6.3): a phased sample `0|1` under G and `0|1` under T at one
  position must not yield a set holding both (nor `1|0` + `1|0`, nor `1|1` + `0|1`, nor three
  ALTs on one slot); it is a carrier of {G} and of {T}, neither CONFIRMED, while a sample
  `1|0` under G and `0|1` under T keeps {G} and {T} CONFIRMED. Add a second column carried
  on the conflicting slot, so the sets with it are checked too.
* The observed PUTATIVE fix (6.3): an unphased `1/2` carrier yields {A}, {C}, never {A, C};
  the same above `max_putative` (the union fallback).
* The population summary of a k ≥ 2 row through a slot holding two ALTs of a row position
  (6.3), and the SNP+SNP companion on a row with two ALTs of one position (6.5): no
  co-occurrence.
* End to end, the run of 6.3 as an acceptance test: on a 2.10.0 index, `complete-search
  --per-sample` with guide `TTATTCCAGACAGATTATGC` (NRG) at chr22:10,661,610 must report no
  row naming two ALTs of one position, in the per-haplotype report, the alternative
  alignments and the companions, and list NA18869 under the G and the T rows, neither as
  CONFIRMED.
  The kit — hg38 chr22:1-10,700,000 and the 4,433 HGDP + 1000 Genomes records up to
  10,699,000 — builds in 15 s and searches in 5 min.
* SNP+indel co-occurrence on a column with two ALTs where the guide uses the non-first one
  (6.4), dictionary and Tier-1 sources.
* An end-to-end guide at a real split site whose earlier record carries a common ALT, e.g.
  chr20:12296886 (G>A, 1,451 copies, absent from the 2.8.x Tier-1 store; G>C, 2 copies,
  present), run through `complete-search` before and after the rebuild.
* A guard test for 6.1: tiers and index from different enrichments are refused.
* `Q` (6.3, 6.4): a target through a `Q` column, alone and beside another variant column,
  yields each of the position's three ALTs in both SNP paths, and the overlay writes `Q` on
  a fake-indel flank base at a four-allele position.

## 7. Decisions of the owner

Taken for CRISPRitz 2.10.0 and CRISPRme 2.9.0 (CRISPRme plan, section 3, B1-B4):

1. **Four-base positions: the fifth symbol `Q`** (B1), so all four alleles are searchable
   (section 2.2). The branch kept the reference and the two most frequent ALTs (by AF),
   leaving out the rarest; at 2,448 chr22 positions with a rarest ALT of AC ≤ 59, the
   measured stake was small, but `Q` loses nothing and costs an index format flag
   (docs/INDEX_FORMAT.md section 4.1), not a rebuild beyond the one already needed.
2. **The same ALT on two records: merge the carriers** (B2; section 2.3), where the branch
   kept the later record's entry (2.9.1's rule, narrowed to the ALT). It happens: 0 such
   positions on chr22 of HGDP + 1000 Genomes, 5 on chr22 of HPRC after `bcftools norm -m
   -any` (4 holding only that ALT, e.g. 23,857,075 T and T; 1 with a third record,
   38,152,965 A, A and C). Measured for the phase: with the later entry kept, HPRC chr22 had
   no slot holding two ALTs (section 4); merging the two A records' carriers (a slot carries
   A if either record puts it there) creates 37 sample-positions with two ALTs on one slot,
   all at 38,152,965, 8 of them with more ALT copies than slots — e.g. HG02071, `1|1` under
   the first A record, `0|0` under the second, `1|1` under C, whom the later entry left a
   carrier of C only. Such conflicts follow decision 4.
3. **Records listing several ALTs are left to CRISPRme** (B3): their entries keep their
   multi-index genotypes (`1|2`), which the registry's `alt_index="1"` and the legacy
   haplotype test mishandle (6.2, 6.3); every CRISPRme reader learns to handle them (store
   the ALT index per allele) rather than CRISPRitz normalising them, so its output for such
   inputs stays stable. Only an entry merged from two records (decision 2) is in split form.
4. **Two ALTs on one haplotype slot are PUTATIVE** (B4), resolved in CRISPRme (6.3), which
   has the window and every ALT's carriers in hand: the enricher passes them through
   verbatim (independently phased split records: 672 sample-positions on HGDP + 1000
   Genomes chr22, none on HPRC chr22), and merging an allele's carriers makes some (decision
   2: 37 on HPRC chr22). Rejected: writing such a sample's genotypes unphased in the
   enricher (it rewrites genotypes the VCF states and still needs the PUTATIVE fix of 6.3),
   and re-phasing each multi-allelic site jointly upstream (a re-phased panel). On the
   CRISPRme side the phase label is per haplotype set, so a set a downgraded sample shares
   turns PUTATIVE for all its carriers (6.3: {G} at chr22:10,661,629, 2,126 carriers,
   through one sample); a phase state per carrier would keep theirs CONFIRMED, at the cost of
   a report format change.

## 8. Published indexes and stores to rebuild

Every variant index built from a panel VCF holding more than one SNP record at a position
(or a record listing three SNP ALTs, section 2.6), and everything derived from its
enrichment. Reference-only indexes (`NRG_3_hg38`, `NRR_3_hg38`, `NGG_3_hg38`, `NNN_3_hg38`,
`NNNNNNNN_3_hg38`) are not affected by the enricher, but are rebuilt too, in index format 2
(docs/INDEX_FORMAT.md; owner decision C1: one format on the Hub, and the sites at contig
ends and N runs format 2 adds).

Panels, each measured on its published chr22 source VCF (`vcfs/<panel>/` of
`lucapinello/crisprme-data`): **HGDP + 1000 Genomes** (`hg38_HGDP1kGP`) — 43,389 positions
with two or more SNP records. **mega** (`hg38_mega`) — 16,103 PASS positions with two or more
SNP ALTs (1,370 with four bases). **HPRC** (`hg38_HPRC`) — the vcfwave VCF holds 4,470
positions with two or more SNP ALTs (61 with four bases), and `hprc_build.sh` splits them
with `bcftools norm -m -any` before enrichment; an HPRC registry compiled from a pre-2.10.0
dictionary (ml008 `allidx/`) holds 535,488 chr22 positions, none with two ALTs. All three
are affected.

`lucapinello/crisprme-data`, `indexes/`:

| Archive | Rebuild | Members that change |
|---|---|---|
| `NRG_3_hg38+hg38_HGDP1kGP` | yes | index bins; `_INDELS` bins (overlay); `Dictionaries/registry_hg38_HGDP1kGP`; `manifest.json` |
| `genotypes_hg38_HGDP1kGP` | yes | the Tier-1 store |
| `NRR_3_hg38+hg38_HGDP1kGP` | yes | as NRG |
| `NGG_3_hg38+hg38_HGDP1kGP` | yes | as NRG |
| `NNN_3_hg38+hg38_HGDP1kGP` | yes | as NRG |
| `NRG_3_hg38+hg38_mega` | yes | index bins; `_INDELS` bins; `manifest.json` (registry unchanged) |
| `NRR_3_hg38+hg38_mega` | yes | as NRG mega |
| `NRG_3_hg38+hg38_HPRC`, `NRR_3_hg38+hg38_HPRC` | yes | as NRG HGDP1kGP, and `Dictionaries/log_indels_hg38_HPRC` (its chrY log, below) |
| `genotypes_hg38_HPRC` | yes | the Tier-1 store |
| `SHA256SUMS` | yes | every replaced archive |

`lucapinello/crisprme-data-general-index`, `indexes/`: `NNNNNNNN_3_hg38+hg38_HGDP1kGP` and
its `genotypes_hg38_HGDP1kGP` (the same companion), and `SHA256SUMS`.

`lucapinello/crisprme-data`, `build-inputs/`: `hg38_HGDP1kGP.tar.gz` and `hg38_mega.tar.gz`
(the enriched genome `Genomes/hg38+hg38_<ds>/` and the overlaid fake-indel genome
`Genomes/hg38+hg38_<ds>_INDELS/`) and their `.SHA256SUMS`.

Unchanged by the enricher (measured on chr22 or derived from the VCF directly):
`indel_genotypes_<vcf>/`, `indel_af_hg38_mega/`, `samplesIDs/`, the source VCFs under
`vcfs/`, and the `log_indels_<vcf>/` logs of HGDP + 1000 Genomes and mega. Not HPRC's: a
log's AF column changes wherever a record lists `INFO/AF` at another place than its VCF's
first record, as each record's AF is now read by its key (section 2.5), and of the 24 HPRC
VCFs (the normalised build input) chrY's has such records (96,853 of its 108,023 indel
records), the 23 others none. With 2.10.0, `log_indels_hg38_HPRC/logchrY` gains the AF of
91,950 of its 102,665 rows (the published one leaves exactly those empty), every other
column identical and each AF its record's (0 mismatches against the VCF), and 47,911 of the
chrY dictionary's 109,033 values change in their AF alone. None of the 23 HGDP + 1000
Genomes VCFs has such a record, and mega's VCF lists no `AF` key (chr22: only `AF_<panel>`
and `AF_max`), so their logs keep their AF; on chr22 the logs of all three panels are
byte-identical to 2.9.1's. Rebuild the logs with the rest rather than carry the published
ones over. A classic (non-dict-less) publish's
`Dictionaries/dictionaries_<vcf>/` changes and must be rebuilt with the rest.

Then: `docs/INDEX_PROVENANCE.md` and `docs/PRECOMPUTED_INDEXES.md` record the rebuild and the
CRISPRitz version, and the release notes say a 2.10.0-built index must not be mixed with
2.8.x tiers (6.1).

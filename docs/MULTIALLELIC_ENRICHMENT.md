# Multi-allelic sites in `add-variants` (CRISPRitz 2.9.2)

What changed in the enricher, what it measured on a real panel, and what CRISPRme must
change, rebuild and test to use it. Written on branch `feat/multiallelic-enrich` (based
on `v2.9.1`). The CRISPRme code cited is the tree at `a5379e0` of the CRISPRme+ 2.8.1
integration branch; line numbers refer to that tree.

**Decided by the owner:** these changes, and the index rebuild they require, go into
CRISPRitz 2.9.2 and the CRISPRme release that pins it, **not** into CRISPRme+ 2.8.1.

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

**From 2.9.2 every ALT is enriched**, in any record order, with one exception that the
search engine forces (section 2.2): at a position whose reference base and ALTs cover all
four bases, the enriched genome leaves out the rarest ALT. The dictionary lists every ALT.

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

### 2.2 Four bases: the rarest ALT is left out of the code

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
| 2.9.2 | either | `V` (A, C, G; T, AF 0.001, left out) | 1, 0 mismatches |

So at such a position the code holds the reference base and its two most frequent ALTs,
leaving out the rarest. An ALT's frequency is the highest of its `INFO/AF` and every
`INFO/AF_<x>` value (`_rank_af` / `rankAf`): a genotyped panel carries `AF`; a sites-only
merge carries only per-dataset fields — CRISPRme's mega has `AF_<dataset>` and `AF_max` and
no `AF`, and ranking by `AF` alone would rank every ALT 0 (on mega's chr22, 1,370 PASS
positions hold four bases, and at 856 of them base order would have left out a more frequent
ALT than the rarest by `AF_max`, at 599 one with `AF_max` ≥ 0.01). A tie leaves out the later
of A, C, G, T; an ALT seen on several records counts with its highest value; a value that is
missing or not a plain decimal number reads as 0. This ranks only: the dictionary's AF field
is still `INFO/AF`. The decision needs every ALT's value, so it is taken once every record is
read (`resolveFourBasePositions`); each single-base ALT's value is kept in three compact
arrays meanwhile (~17 bytes per ALT). It does not depend on the dictionary, so mode `no`
(genome only) codes exactly what mode `yes` codes (unit-tested). One line on **stdout**
counts these positions per chromosome (stderr is fatal to a CRISPRme stage).

This also fixes a latent defect of 2.9.1: a record listing three SNP ALTs wrote `N`.

On chr22 of HGDP + 1000 Genomes, 2,448 positions carry three SNP ALTs. The rarest ALT has
AC = 1 at 1,193 of them, AC ≥ 10 at 110, AC ≥ 41 (0.5 % of 8,182 haplotypes) at 3, and at
most 59; the two rarest tie at 623 (section 4). An engine that could store a four-base
code would lift the exception (section 7, decision 1).

### 2.3 The SNP dictionary (`my_dict_<chr>.json`)

The value of `"<chr>,<pos>"` lists **every ALT** of the position: one entry per ALT joined
by `$`, the form a record listing several ALTs has always produced, so every reader of
multi-ALT values reads it unchanged. Each entry is `<sample:GT,...>;<REF>,<ALT>;<rsID>;<AF>`
built from **its own record**, exactly as 2.9.1 built a lone record's value: its carriers
(sorted), their genotypes verbatim (for a split record `1` is that entry's ALT, which is
`tier0_registry`'s default `alt_index`), its rsID (`;` normalised to `,`) and its own AF.
Entries are in file order. A second record for an ALT already listed (the same allele twice)
replaces that ALT's entry in place, as the whole value used to be replaced; CRISPRme's
Tier-1 compiler refuses a repeated `(pos, alt)` (`tier1_genotypes.py:283-287`), so the
dictionary never holds one. Keys and key order are unchanged. A key holding one record is
byte-identical to 2.9.1.

### 2.4 The indel log and the fake-indel genome

Unchanged. `indel_to_fasta` cuts every contig from the plain reference (`genomeStr`), never
from the enriched genome, so a merged SNP code next to or at an indel does not reach a
contig, its log row or its FAKEPOS (unit test `test_indel_outputs_unmoved_by_split_snps`;
`verify_enricher.sh` cases 19, 21, 28; byte-identical on chr22 of all three panels, section
4). The multi-allelic indel fix of 2.9.1 (each allele's contig from the reference window) is
untouched and still covered by case 7. CRISPRme's overlay (`overlay_indel_snps.py`) copies
the new codes onto the contigs' flanks afterwards (section 6.4).

### 2.5 FILTER and AF act per ALT

A record whose FILTER is neither `PASS` nor `.` contributes nothing, to the code or the
dictionary, and the PASS records at its position are still merged (test
`test_non_pass_record_contributes_nothing`, case 23). Each ALT keeps its own AF in the
dictionary and in the four-base rule. A minor-allele-frequency filter applied before
enrichment acts per ALT only on split records — measured with bcftools 1.24 on a site
`G>A,C,T` with AF 0.3, 0.0005, 0.02 and `bcftools view -i 'INFO/AF>0.001'`:

| input to `add-variants` | 2.9.1 code / dict ALTs | 2.9.2 code / dict ALTs |
|---|---|---|
| split, unfiltered | `K` / T | `D` (C left out as rarest) / A, C, T |
| split, then filtered (C dropped) | `K` / T | `D` / A, T |
| one record, filtered (kept whole: one ALT passes) | **`N`** / A, C, T | `D` / A, C, T |

So filter **after** splitting, as the published panels are built.

### 2.6 Inputs without split multi-allelic SNP sites

Byte-identical to 2.9.1 in all four outputs and the run log: measured on the 17-case
`verify_enricher.sh` fixture of 2.9.1 with its two split sites reduced to one record each
(both ports of 2.9.1 and of this branch: identical md5 for every output). The only other
changes of output are a record listing three SNP ALTs (2.2) and a SNP outside the
chromosome (2.7).

### 2.7 A SNP outside the chromosome

A SNP record whose POS is 0 or past the last base (a VCF of another assembly) now stops the
enricher with `ERROR: <chr>:<pos> is outside the chromosome (<n> bp)` on stderr and exit
status 1, in both ports; `add-variants` then fails naming the chromosome. In 2.9.1 the
Python port raised an `IndexError` past the end and coded the **last** base for POS 0, and
the C++ port read and wrote out of bounds; the new per-position mask would have made that
write unconditional.

## 3. Cost

One byte per base for the union and ~17 bytes per single-base ALT for the four-base rule,
in both ports; the merge itself runs only at a key already present. Measured on chr22
(section 4).

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
(the index-only install's, compiled from a pre-2.9.2 dictionary) holds 1,711,938 positions on
chr22 — exactly the positions with a SNP record — none with two ALTs.

Before vs after (compiled ports; the Python port, below, is byte-identical to the compiled
one):

| output | result |
|---|---|
| enriched genome | **43,389 bases change, exactly the multi-SNP positions**; every base equals the rule: 2.9.1 = reference + last record's ALT, 2.9.2 = the merged rule of 2.1-2.2, at all 1,770,802 positions with a REF-length-1 record (0 mismatches each); no base outside SNP positions differs from the reference in either; `N` written for a variant: 0 and 0 |
| codes | all 2-base → 3-base: Y→H 9,349, S→V 8,623, K→D 8,262, Y→B 6,943, K→B 3,838, W→D 2,175, R→V 2,079, W→H 876, K→V 526, W→V 317, K→H 210, Y→V 191 |
| non-last ALTs | all 45,837 absent from the 2.9.1 code; 44,633 present in the new code; 1,204 are the rarest ALT of a four-base position |
| four-base positions | 2,448 (the stdout note says so). At 1,244 the rarest ALT is the last record, which 2.9.1 coded and 2.9.2 leaves out for the two others (AC = 1: 626, 2-9: 571, 10-40: 46, ≥ 41: 1); at 1,204 it is an earlier record (AC = 1: 567, 2-9: 574, 10-40: 61, ≥ 41: 2) |
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

## 5. Tests on the branch

* `test/test_enricher_multiallelic.py` — 21 unit tests, each run against **both** ports
  (42 in all), with exact expected outputs: two and three ALTs in every order, split ==
  joined, the four-base rule (AF order, ties, AF spellings, ranking by `AF_<x>`/`AF_max`
  without `AF`, an ALT on two records, a record listing three ALTs), reference `N` and
  IUPAC bases, mode `no` == mode `yes`, the dictionary entries (carriers, rsIDs, AFs,
  haploid / unphased / missing / `GT:DP`
  genotypes, records without carriers), a non-PASS record, the same ALT twice, a split
  record beside a multi-ALT record, `*`, non-adjacent records, the indel outputs, a
  position with one record, and a SNP outside the chromosome. Against 2.9.1 38 of the 42
  fail; the 4 that pass (mode `no` == mode `yes`, a position with one record) pin what does
  not change. In CI (`enricher-ports`).
* `verify_enricher.sh` — grows from 17 to 29 cases (cases 18-29: split sites in every
  order, the same sites as one record, the four-base rule, non-adjacent records, genotype
  shapes, SNP + indel at one position and in an indel's window, a split record beside a
  multi-ALT record, a non-PASS record, the same ALT twice, AF spellings, `AF_<x>` fields, a
  reference IUPAC base and `N`, `*`, the chromosome edges); both ports byte-identical, with
  the file set and the run log. Checked against an independent oracle (89 coded positions,
  79 dictionary keys, 0 mismatches).
* The existing scripts pass unchanged on Linux (g++ 14): `run_minimal_test.sh` (NGG, WTN,
  WTN under AddressSanitizer), `run_golden_test.sh`, `run_pam_filter_test.sh` (31
  assertions), `run_flank_n_test.sh`, `run_short_contig_test.sh`,
  `run_add_variants_status_test.sh`.

## 6. What CRISPRme must change

Nothing in CRISPRme's file formats has to change: the Tier-0 registry and the Tier-1 store
already hold one record per `(pos, alt)` and several per position, and every dictionary
reader already splits a value on `$`. What changes is **content** (positions that held one
ALT now hold two or three), so the data must be rebuilt (section 8), and a few places that
assumed one ALT per position, or that were only ever exercised with one, must be fixed.

### 6.1 Pin, provenance, and a guard against mixing

* Pin CRISPRitz `2.9.2`: `Dockerfile:31` (`ARG crispritz_ref=v2.9.1`), `environment.yml`
  and `install_from_source.sh` (which name 2.9.1), the Bioconda `crisprme` recipe.
* The enriched genome, its TST index, the fake-indel genome and its `_INDELS` index, the
  per-sample dictionaries and the Tier-0/Tier-1 stores must all come from **one**
  enrichment. Mixed, the code and the genotype source disagree at ~2.5 % of SNP positions:
  an old dictionary under a new genome has no carriers for a coded ALT (rows needing it are
  not emitted), a new dictionary under an old genome offers ALTs the search never matched. Stamp the enrichment in the tier manifests
  (`registry_<vcf>/reg_<chrom>.idx`, whose manifest `tier0_registry._write_registry`
  already stamps with `data_type`/`phased`, line 1158, and `genotypes_<vcf>/gt_<chrom>.idx`,
  written by `tier1_genotypes.compile_genotypes`) and in the index's `manifest.json`, e.g.
  `"enricher": {"crispritz": "2.9.2", "multiallelic": "merged"}`, and refuse (or warn on)
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
  entries (a slot holding the other ALT reads `0`). No change; but note it is **wrong for
  entries of a record listing several ALTs** (`1|2` puts the sample on both slots of both
  ALTs) — pre-existing, not reachable on the published panels; the observed path below
  resolves it per record and could be reused here. This path (and its population-level
  greedy, `by_pos`, from line 1027) substitutes **every** ALT the dictionary lists, without
  checking the target's IUPAC column, so at a four-base position it also tries the rarest
  ALT the code left out; such a row passes only if it fits the mismatch budget on its own.
  Filtering by the column, as the observed path does, would make the two paths agree.
* Observed-haplotype path, `_iupac_decomposition_observed` (lines 607-): the alt index is
  resolved per record from its carriers (line 675, `observed_haplotypes._alt_index_for_record`)
  — `1` for split-record entries — and the target's own IUPAC column filters the ALTs (#139,
  line 655), which also drops the rarest ALT at a four-base position, consistently with
  the search. CONFIRMED (phased) sets are right: measured with CRISPRme's own module, a
  phased carrier `1|0` under A and `0|1` under C gives the two sets {A} and {C}.
* **Must change — `observed_haplotypes._per_sample_variant_sets`, PUTATIVE branch (lines
  223-249).** An unphased carrier of two ALTs of one position (a `1/2` heterozygote, split
  into `1/0` under one ALT and `0/1` under the other) gets every subset of its carried
  variants as a candidate, **including the set holding both ALTs of that column**: measured,
  the enumerator emits `{(3,A),(3,C)}` with a sequence holding C and info listing both
  variants — a row naming two alleles at one position. Unreachable while the dictionary held
  one ALT per position; reachable now on every **unphased** panel. Fix: build the PUTATIVE
  combinations with at most one ALT per column (`itertools.product` over each carried
  column's {reference, ALTs} minus the empty choice), or drop combinations whose columns
  repeat. Phased panels (HGDP + 1000 Genomes, HPRC) do not reach it.
* `population_summary` (`_phased_decide`/`_unphased_decide`, lines 648-698, `alt_index="1"`)
  and `twopass_emit` (one candidate chosen per column): no change.
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
  2.9.2 pin. Left over there: `forward_snp_tokens` still picks `alts[0]` but has no caller —
  delete it. The column test `if col.ref not in col.code` stays true: the four-base rule
  always keeps the reference base.
* `overlay_indel_snps.py`: copies each flank base's code from the enriched genome
  (`_IUPAC_AMBIG`, line 34, includes `B`/`D`/`H`/`V`; the enricher never writes `N` for a
  variant). No change; the overlaid fake-indel genome and its `_INDELS` index must be
  **rebuilt** (section 4: measured which flank bases change).
* `build_indel_genotypes.py` (keyed `<pos>_<REF>_<ALT>`) and the log-driven carriers: not
  affected.

### 6.5 The report

No format change: a row's variant columns list the ALT each haplotype holds, and two ALTs of
one position appear in two rows, as they would for a record listing several ALTs
(`resultIntegrator.py` splits the tokens per row, lines 455-520). Expect more variant rows
and more carriers at the affected positions; the changelog and `METHODS.md` must say that
earlier releases reported only the last ALT of a split site, and the four-base rule.

### 6.6 Tests CRISPRme should add

* A Tier-0/Tier-1 compile of a merged value (two ALTs, carriers on opposite slots): per-ALT
  AC/AN and carriers (the measurement in 6.2, as a unit test).
* The observed PUTATIVE fix (6.3): an unphased `1/2` carrier yields {A}, {C}, never {A, C}.
* SNP+indel co-occurrence on a column with two ALTs where the guide uses the non-first one
  (6.4), dictionary and Tier-1 sources.
* An end-to-end guide at a real split site whose earlier record carries a common ALT, e.g.
  chr20:12296886 (G>A, 1,451 copies, absent from the 2.8.x Tier-1 store; G>C, 2 copies,
  present), run through `complete-search` before and after the rebuild.
* A guard test for 6.1: tiers and index from different enrichments are refused.

## 7. Open decisions for the owner

1. **Four-base positions.** This branch keeps the reference and the two most frequent ALTs
   (by AF). The alternative is an engine change: a fifth symbol that `buildTST`'s writer,
   the loader and `searchTST`'s comparison read as `1111`, so `N` stays a gap — an index
   format change. At 2,448 chr22 positions with a rarest ALT of AC ≤ 59, the measured stake
   is small.
2. **The same ALT on two records** keeps the later record's entry (2.9.1's rule, narrowed
   to the ALT), so at a position holding only that ALT the output is 2.9.1's. Merging the two
   records' carriers is the other choice. It happens: 0 such positions on chr22 of HGDP +
   1000 Genomes, 5 on chr22 of HPRC after `bcftools norm -m -any` (4 holding only that ALT,
   e.g. 23,857,075 T and T; 1 with a third record, 38,152,965 A, A and C).
3. **Records listing several ALTs** keep their multi-index genotypes, which the registry's
   `alt_index="1"` and the legacy haplotype test mishandle (6.2, 6.3). Normalising them to per-ALT genotypes in the enricher would
   change the output for such inputs; leaving them to CRISPRme (as `r281f-indel` does for
   the indel path) keeps CRISPRitz's output stable.

## 8. Published indexes and stores to rebuild

Every variant index built from a panel VCF holding more than one SNP record at a position,
and everything derived from its enrichment. Reference-only indexes are not affected:
`NRG_3_hg38`, `NRR_3_hg38`, `NGG_3_hg38`, `NNN_3_hg38`, `NNNNNNNN_3_hg38`.

Panels, each measured on its published chr22 source VCF (`vcfs/<panel>/` of
`lucapinello/crisprme-data`): **HGDP + 1000 Genomes** (`hg38_HGDP1kGP`) — 43,389 positions
with two or more SNP records. **mega** (`hg38_mega`) — 16,103 PASS positions with two or more
SNP ALTs (1,370 with four bases). **HPRC** (`hg38_HPRC`) — the vcfwave VCF holds 4,470
positions with two or more SNP ALTs (61 with four bases), and `hprc_build.sh` splits them
with `bcftools norm -m -any` before enrichment; an HPRC registry compiled from a pre-2.9.2
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
| `NRG_3_hg38+hg38_HPRC`, `NRR_3_hg38+hg38_HPRC` | yes | as NRG HGDP1kGP |
| `genotypes_hg38_HPRC` | yes | the Tier-1 store |
| `SHA256SUMS` | yes | every replaced archive |

`lucapinello/crisprme-data-general-index`, `indexes/`: `NNNNNNNN_3_hg38+hg38_HGDP1kGP` and
its `genotypes_hg38_HGDP1kGP` (the same companion), and `SHA256SUMS`.

`lucapinello/crisprme-data`, `build-inputs/`: `hg38_HGDP1kGP.tar.gz` and `hg38_mega.tar.gz`
(the enriched genome `Genomes/hg38+hg38_<ds>/` and the overlaid fake-indel genome
`Genomes/hg38+hg38_<ds>_INDELS/`) and their `.SHA256SUMS`.

Unchanged by the enricher (measured on chr22 or derived from the VCF directly): the
`log_indels_<vcf>/` logs, `indel_genotypes_<vcf>/`, `indel_af_hg38_mega/`, `samplesIDs/`,
the source VCFs under `vcfs/`. A classic (non-dict-less) publish's
`Dictionaries/dictionaries_<vcf>/` changes and must be rebuilt with the rest.

Then: `docs/INDEX_PROVENANCE.md` and `docs/PRECOMPUTED_INDEXES.md` record the rebuild and the
CRISPRitz version, and the release notes say a 2.9.2-built index must not be mixed with
2.8.x tiers (6.1).

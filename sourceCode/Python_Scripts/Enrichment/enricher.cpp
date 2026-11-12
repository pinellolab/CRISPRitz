// C++ port of enricher.py (CRISPRitz variant-enrichment / "add-variants" step).
//
// Behaviorally byte-identical to sourceCode/Python_Scripts/Enrichment/enricher.py.
// Reproduces the same four outputs (enriched FASTA, my_dict_<chr>.json,
// fake<chr>.fa, log<chr>.txt) including the original script's cwd side effects.
//
// Build: $(CXX) -std=c++11 -O3 -fopenmp ... enricher.cpp -o enricher -lz
//
// Usage: enricher <vcf.gz> <chrom.fa> <argv3prefix> <yes|no> <argv5path>

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <set>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include <zlib.h>

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------

static std::vector<std::string> splitStr(const std::string &s, char delim) {
    std::vector<std::string> out;
    size_t start = 0;
    while (true) {
        size_t pos = s.find(delim, start);
        if (pos == std::string::npos) {
            out.push_back(s.substr(start));
            break;
        }
        out.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return out;
}

// ALT allele indices carried by a VCF GT field.
// Splits on '|' and '/', so phased, unphased, haploid ("1") and
// partially-missing (".|1") genotypes are all handled. Reference (0) and
// missing ('.') alleles are skipped, so the result holds only ALT indices.
// Matching on these integers -- instead of searching the raw GT text for an
// allele number -- is what keeps "0|12" from matching allele 1.
static std::set<int> gtAltAlleles(const std::string &gt) {
    std::set<int> out;
    std::vector<std::string> parts = splitStr(gt, '|');
    for (size_t i = 0; i < parts.size(); ++i) {
        std::vector<std::string> sub = splitStr(parts[i], '/');
        for (size_t j = 0; j < sub.size(); ++j) {
            const std::string &a = sub[j];
            if (a.empty()) continue;
            bool numeric = true;
            for (size_t c = 0; c < a.size(); ++c) {
                if (a[c] < '0' || a[c] > '9') { numeric = false; break; }
            }
            if (!numeric) continue;  // '.' (missing) or anything unexpected
            int v = atoi(a.c_str());
            if (v > 0) out.insert(v);  // 0 is the reference allele
        }
    }
    return out;
}

// Python str.strip(): removes leading/trailing ASCII whitespace.
static std::string pyStrip(const std::string &s) {
    size_t b = 0, e = s.size();
    while (b < e && (s[b] == ' ' || s[b] == '\t' || s[b] == '\n' ||
                     s[b] == '\r' || s[b] == '\f' || s[b] == '\v'))
        ++b;
    while (e > b && (s[e - 1] == ' ' || s[e - 1] == '\t' || s[e - 1] == '\n' ||
                     s[e - 1] == '\r' || s[e - 1] == '\f' || s[e - 1] == '\v'))
        --e;
    return s.substr(b, e - b);
}

// ---------------------------------------------------------------------------
// gzip line reader
// ---------------------------------------------------------------------------

struct GzReader {
    gzFile fp;
    std::string buf;
    size_t pos;
    bool eof;
    GzReader(const char *path) : pos(0), eof(false) {
        fp = gzopen(path, "rb");
        if (!fp) {
            fprintf(stderr, "ERROR: cannot open %s\n", path);
            exit(1);
        }
    }
    ~GzReader() {
        if (fp) gzclose(fp);
    }
    // Returns false at EOF with no more data. Line includes no trailing '\n'
    // (we strip it here since callers always strip/split anyway, except the
    // VCF header read which in python keeps the split-on-tab of a raw line).
    bool getRawLine(std::string &line) {
        line.clear();
        while (true) {
            if (pos >= buf.size()) {
                if (eof) return !line.empty();
                char tmp[1 << 16];
                int n = gzread(fp, tmp, sizeof(tmp));
                if (n <= 0) {
                    eof = true;
                    if (line.empty()) return false;
                    return true;
                }
                buf.assign(tmp, tmp + n);
                pos = 0;
            }
            while (pos < buf.size()) {
                char c = buf[pos++];
                if (c == '\n') return true;
                line.push_back(c);
            }
        }
    }
};

// ---------------------------------------------------------------------------
// IUPAC mapping: base SET -> ambiguity letter (order-independent).
// Canonicalize the set of {A,C,G,T} present into a 4-bit mask.
// ---------------------------------------------------------------------------

static inline int baseBit(char c) {
    switch (c) {
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 4;
        case 'T': return 8;
        default:  return 0;  // non-ACGT contributes nothing recognizable
    }
}

// Whether a reference byte is a letter a SNP is coded over: A/C/G/T and the IUPAC
// ambiguity codes. N is a gap, and so is any other byte (buildTST and searchBruteForce
// read it as N); Q is refused when the FASTA is read. enricher.py's _REFERENCE_CODES.
static inline bool isReferenceCode(char c) {
    return c != '\0' && strchr("ACGTRYSWKMBDHV", c) != NULL;
}

// Maps a mask of ACGT bits to the IUPAC letter reproduced by enricher.py.
static char iupacFromMask(int mask) {
    switch (mask) {
        // A SINGLE base is NOT a key in python's iupac_code: its values are the
        // 2-, 3- and 4-base combinations only (12 + 24 + 24 entries, none of
        // length 1), so when the allele set collapses to one base python's lookup
        // finds nothing and leaves the genome UNTOUCHED. Returning the base here
        // made the caller WRITE it, clobbering an IUPAC code that a SNP record at
        // the same position had already placed.
        //
        // The set collapses to one base on every indel record: REF is length 1, so
        // SNPsProcess does not bail out, and the longer ALT is skipped by the
        // `size() > 1` test, leaving only the reference base. A position carrying
        // BOTH a SNP and an indel -- as separate biallelic records, which is how
        // every CRISPRme-supplied VCF is normalised -- therefore lost its SNP.
        // Measured on a real 120,000-record chr22 slice: 176 such positions were
        // the ONLY bytes differing from enricher.py in 50.8 MB.
        // verify_enricher.sh case 8 pins it.
        case 1:  // A
        case 2:  // C
        case 4:  // G
        case 8:  // T
            return 0;
        case (1 | 4):     return 'R';  // A,G
        case (2 | 8):     return 'Y';  // C,T
        case (4 | 2):     return 'S';  // G,C
        case (1 | 8):     return 'W';  // A,T
        case (4 | 8):     return 'K';  // G,T
        case (1 | 2):     return 'M';  // A,C
        case (2 | 4 | 8): return 'B';  // C,G,T
        case (1 | 4 | 8): return 'D';  // A,G,T
        case (1 | 2 | 8): return 'H';  // A,C,T
        case (1 | 2 | 4): return 'V';  // A,C,G
        // A,C,G,T: Q, which the search reads as any base -- not N, an assembly gap that
        // matches nothing (enricher.py's _IUPAC_BY_MASK[15])
        case (1 | 2 | 4 | 8): return 'Q';
        default: return 0;  // no matching key -> genome unchanged (python: no key hit)
    }
}

// ---------------------------------------------------------------------------
// case-insensitive replacement of the FIRST occurrence of `pat` in `s`,
// reproducing re.sub(pat, repl, s, count=1, flags=IGNORECASE).
// The VCF ref/alt strings here are plain literals (no regex metachars in
// nucleotide sequences), so a literal case-insensitive find is equivalent.
// ---------------------------------------------------------------------------

static std::string reSubFirstCI(const std::string &pat, const std::string &repl,
                                const std::string &s) {
    if (pat.empty()) {
        // re.sub with empty pattern inserts at position 0; not expected for
        // nucleotide refs (always length >= 1). Guard anyway.
        return repl + s;
    }
    size_t n = s.size(), m = pat.size();
    for (size_t i = 0; i + m <= n; ++i) {
        bool match = true;
        for (size_t j = 0; j < m; ++j) {
            char a = s[i + j], b = pat[j];
            if (a >= 'a' && a <= 'z') a -= 32;
            if (b >= 'a' && b <= 'z') b -= 32;
            if (a != b) { match = false; break; }
        }
        if (match) {
            return s.substr(0, i) + repl + s.substr(i + m);
        }
    }
    return s;  // no match -> unchanged, matching re.sub behavior
}

// ---------------------------------------------------------------------------
// globals mirroring the python module-level state
// ---------------------------------------------------------------------------

static std::string genomeStr;
static std::string genomeHeader;    // includes leading '>' and trailing '\n'
static std::string currentChr;
static std::vector<std::string> VCFheader;

// enriched genome is edited in place on a copy of genomeStr
static std::string genomeList;

// The A/C/G/T bases coded so far at each position, as a 4-bit mask (0 = no SNP
// record here yet) -- enricher.py's snp_mask. A split multi-allelic site arrives as
// one record per ALT, so each record adds its ALT to what the earlier ones coded.
// Kept beside genomeList because a letter there cannot tell a code written here from
// an ambiguity code the reference FASTA already holds.
static std::vector<unsigned char> snpMask;
// positions coded Q (all four bases), counted for the stdout note
static long fourBaseCount = 0;

// SNP dict: preserve insertion order (python dict) -> use vector of pairs,
// with a lookup map for "already present" merge semantics.
static std::vector<std::pair<std::string, std::string>> chr_dict_snps;
static std::map<std::string, size_t> chr_dict_index;
// (key, "REF,ALT") -> the allele index of that entry's ALT in its record, for an entry
// of a record listing several ALTs (1, a split record's, is not stored) --
// enricher.py's entry_allele
static std::map<std::pair<std::string, std::string>, int> entryAllele;

static std::string joinStr(const std::vector<std::string> &parts, char sep) {
    std::string out;
    for (size_t k = 0; k < parts.size(); ++k) {
        if (k) out += sep;
        out += parts[k];
    }
    return out;
}

// The four ';' fields of one dict entry "<samples>;<REF>,<ALT>;<rsID>;<AF>", split from
// the right like enricher.py's entry.rsplit(';', 3), so a ';' in a sample name cannot
// move them. Every entry this program builds holds at least three ';'.
static std::vector<std::string> entryFields(const std::string &e) {
    size_t p3 = e.rfind(';');
    size_t p2 = (p3 == std::string::npos || p3 == 0) ? std::string::npos : e.rfind(';', p3 - 1);
    size_t p1 = (p2 == std::string::npos || p2 == 0) ? std::string::npos : e.rfind(';', p2 - 1);
    std::vector<std::string> f(4);
    if (p1 == std::string::npos) {  // unreachable for entries built here
        f[1] = e;
        return f;
    }
    f[0] = e.substr(0, p1);
    f[1] = e.substr(p1 + 1, p2 - p1 - 1);
    f[2] = e.substr(p2 + 1, p3 - p2 - 1);
    f[3] = e.substr(p3 + 1);
    return f;
}

static bool allDigits(const std::string &a) {
    if (a.empty()) return false;
    for (size_t c = 0; c < a.size(); ++c)
        if (a[c] < '0' || a[c] > '9') return false;
    return true;
}

// A GT split into alleles and the separators between them, as enricher.py's
// re.split(r'([|/])', gt): even entries alleles, odd entries '|' or '/'.
static std::vector<std::string> gtTokens(const std::string &gt) {
    std::vector<std::string> t(1);
    for (size_t c = 0; c < gt.size(); ++c) {
        if (gt[c] == '|' || gt[c] == '/') {
            t.push_back(std::string(1, gt[c]));
            t.push_back("");
        } else {
            t.back().push_back(gt[c]);
        }
    }
    return t;
}

// GT gt with allele k written 1 and every other called allele 0 -- enricher.py's
// _gt_for_allele: the form of a split record's genotypes, where 1 is its one ALT.
static std::string gtForAllele(const std::string &gt, int k) {
    std::vector<std::string> t = gtTokens(gt);
    std::string out;
    for (size_t i = 0; i < t.size(); ++i) {
        if (allDigits(t[i]))
            out += (atol(t[i].c_str()) == k) ? "1" : "0";
        else
            out += t[i];
    }
    return out;
}

// One sample's genotypes of one allele on two records, the earlier first, both in split
// form -- enricher.py's _gt_union: phased alike (as many alleles, every separator '|'), a
// slot holds the allele if either record puts it there, else a called 0 if either has
// one; otherwise the genotype with more copies of the allele (on a tie, the earlier).
static std::string gtUnion(const std::string &gt1, const std::string &gt2) {
    std::vector<std::string> t1 = gtTokens(gt1), t2 = gtTokens(gt2);
    bool alike = t1.size() == t2.size();
    for (size_t i = 0; alike && i < t1.size(); ++i)
        if (t1[i] == "/" || t2[i] == "/") alike = false;
    if (!alike)
        return std::count(t2.begin(), t2.end(), std::string("1")) >
                       std::count(t1.begin(), t1.end(), std::string("1"))
                   ? gt2
                   : gt1;
    std::string out;
    for (size_t i = 0; i < t1.size(); ++i) {
        const std::string &x = t1[i], &y = t2[i];
        if (x == "|")
            out += x;
        else if (x == "1" || y == "1")
            out += "1";
        else if (x == "0" || y == "0")
            out += "0";
        else
            out += x;
    }
    return out;
}

static double afRank(const std::string &af);

// One allele's dict entries from two records (the same REF>ALT twice) -- enricher.py's
// _merge_entries: every carrier of either, genotypes in split form, a carrier of both
// with the union of its two; the rsID and AF of the more frequent record (highest
// INFO/AF; on a tie, the earlier).
static std::string mergeEntries(const std::string &oldE, int kOld, const std::string &newE, int kNew) {
    std::vector<std::string> fo = entryFields(oldE), fn = entryFields(newE);
    std::map<std::string, std::string> gts;  // sample -> GT
    const std::string *samples[2] = {&fo[0], &fn[0]};
    const int ks[2] = {kOld, kNew};
    for (int r = 0; r < 2; ++r) {
        if (samples[r]->empty()) continue;
        std::vector<std::string> carriers = splitStr(*samples[r], ',');
        for (size_t c = 0; c < carriers.size(); ++c) {
            size_t colon = carriers[c].rfind(':');
            std::string name = carriers[c].substr(0, colon);
            std::string gt = gtForAllele(carriers[c].substr(colon + 1), ks[r]);
            std::map<std::string, std::string>::iterator it = gts.find(name);
            if (it == gts.end())
                gts[name] = gt;
            else
                it->second = gtUnion(it->second, gt);
        }
    }
    std::vector<std::string> listed;
    for (std::map<std::string, std::string>::const_iterator it = gts.begin(); it != gts.end(); ++it)
        listed.push_back(it->first + ":" + it->second);
    std::sort(listed.begin(), listed.end());
    bool newer = afRank(fn[3]) > afRank(fo[3]);
    return joinStr(listed, ',') + ";" + fo[1] + ";" + (newer ? fn[2] : fo[2]) + ";" + (newer ? fn[3] : fo[3]);
}

// Store one record's per-ALT entries under "chr,pos", MERGED with what earlier
// records stored there -- enricher.py's _dict_put, see it for the rationale. A new
// key gets the entries joined by '$' (as before); at an existing key an allele already
// listed (the same REF>ALT) gets the carriers of both (mergeEntries), any other ALT is
// appended. alleles[i] is the allele index of entries[i]'s ALT in its record.
static void dictPut(const std::string &key, const std::vector<std::string> &entries,
                    const std::vector<int> &alleles) {
    std::map<std::string, size_t>::iterator it = chr_dict_index.find(key);
    if (it == chr_dict_index.end()) {
        chr_dict_index[key] = chr_dict_snps.size();
        chr_dict_snps.push_back(std::make_pair(key, joinStr(entries, '$')));
        for (size_t k = 0; k < entries.size(); ++k)
            if (alleles[k] != 1) entryAllele[std::make_pair(key, entryFields(entries[k])[1])] = alleles[k];
        return;
    }
    std::string &val = chr_dict_snps[it->second].second;  // merge, keep position
    std::vector<std::string> merged = splitStr(val, '$');
    std::vector<std::string> refalts;
    for (size_t k = 0; k < merged.size(); ++k) refalts.push_back(entryFields(merged[k])[1]);
    for (size_t k = 0; k < entries.size(); ++k) {
        std::string ra = entryFields(entries[k])[1];
        size_t j = 0;
        while (j < refalts.size() && refalts[j] != ra) ++j;
        if (j < refalts.size()) {
            std::pair<std::string, std::string> id(key, ra);
            std::map<std::pair<std::string, std::string>, int>::iterator ka = entryAllele.find(id);
            int kOld = 1;
            if (ka != entryAllele.end()) {
                kOld = ka->second;
                entryAllele.erase(ka);
            }
            merged[j] = mergeEntries(merged[j], kOld, entries[k], alleles[k]);
        } else {
            merged.push_back(entries[k]);
            refalts.push_back(ra);
            if (alleles[k] != 1) entryAllele[std::make_pair(key, ra)] = alleles[k];
        }
    }
    val = joinStr(merged, '$');
}

// fake indel fasta pieces and log rows
static std::vector<std::string> list_fasta_indels;
static std::vector<std::vector<std::string>> log_indels;  // 7 cols each

static bool doYes = false;

// ---------------------------------------------------------------------------
// isPlainNumber(s): enricher.py's _AF_NUMBER, the grammar
// [+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?, checked before strtod, which
// on its own would also take "nan", "inf", hex and leading spaces, where Python's
// float() takes other spellings ("1_000", trailing spaces): the check is what makes
// both ports read every AF string alike.
// ---------------------------------------------------------------------------

static inline bool isDigitC(char c) { return c >= '0' && c <= '9'; }

static bool isPlainNumber(const std::string &s) {
    size_t i = 0, n = s.size();
    if (i < n && (s[i] == '+' || s[i] == '-')) ++i;
    size_t intDigits = 0, fracDigits = 0;
    while (i < n && isDigitC(s[i])) { ++i; ++intDigits; }
    if (i < n && s[i] == '.') {
        ++i;
        while (i < n && isDigitC(s[i])) { ++i; ++fracDigits; }
    }
    if (intDigits == 0 && fracDigits == 0) return false;
    if (i < n && (s[i] == 'e' || s[i] == 'E')) {
        size_t j = i + 1, expDigits = 0;
        if (j < n && (s[j] == '+' || s[j] == '-')) ++j;
        while (j < n && isDigitC(s[j])) { ++j; ++expDigits; }
        if (expDigits == 0) return false;
        i = j;
    }
    return i == n;
}

// afRank(af): an AF value string as a number, to tell which of two records of one
// allele is the more frequent (mergeEntries) -- enricher.py's _af_rank: 0 when it is
// missing or not a plain decimal number.
static double afRank(const std::string &af) {
    return isPlainNumber(af) ? strtod(af.c_str(), NULL) : 0.0;
}

// ---------------------------------------------------------------------------
// SNPsProcess(line): mutates genomeList and line[1] -> 0-based.
// ---------------------------------------------------------------------------

static void SNPsProcess(std::vector<std::string> &line) {
    if (line[3].size() > 1) return;  // not a SNP
    long pos0 = strtol(line[1].c_str(), NULL, 10) - 1;
    line[1] = std::to_string(pos0);  // python sets line[1] to 1-based-1
    if (pos0 < 0 || (size_t)pos0 >= genomeStr.size()) {
        // enricher.py fails here too (it raised IndexError past the end and coded the
        // last base for POS 0); reading and writing out of bounds is no answer.
        fprintf(stderr,
                "ERROR: %s:%ld is outside the chromosome (%ld bp): is the VCF on another "
                "assembly?\n",
                currentChr.c_str(), pos0 + 1, (long)genomeStr.size());
        exit(1);
    }

    // An N is an assembly gap and stays one, as does any other byte that is no base or
    // IUPAC code, which buildTST and searchBruteForce read as N (see enricher.py's
    // SNPsProcess). The dict still lists the record.
    if (!isReferenceCode(genomeStr[pos0])) return;

    // build the base SET: what earlier records coded here (else the reference base)
    // + each alt allele of length 1. Starting from the reference base alone on each
    // record kept only the last record's ALT of a split multi-allelic site; see
    // enricher.py's SNPsProcess. All four bases are Q (iupacFromMask).
    int mask = snpMask[pos0] ? snpMask[pos0] : baseBit(genomeStr[pos0]);
    std::vector<std::string> altAlleles = splitStr(pyStrip(line[4]), ',');
    for (size_t k = 0; k < altAlleles.size(); ++k) {
        if (altAlleles[k].size() != 1) continue;
        mask |= baseBit(altAlleles[k][0]);
    }
    if (mask == 15 && snpMask[pos0] != 15) ++fourBaseCount;
    snpMask[pos0] = (unsigned char)mask;
    char letter = iupacFromMask(mask);
    if (letter != 0) {
        genomeList[pos0] = letter;
    }
    // if letter==0, python loop found no matching value -> genome unchanged.
}

// ---------------------------------------------------------------------------
// afValue(line): raw INFO/AF value string of this record, found by its exact key,
// or "" if the record has none. Mirrors enricher.py's _af_value() -- see that
// function's comment: until 2.9.1 the AF field's place was taken from the first
// record only, and a record listing its INFO fields otherwise got "".
// ---------------------------------------------------------------------------

static std::string afValue(const std::vector<std::string> &line) {
    std::vector<std::string> splitted = splitStr(line[7], ';');
    for (size_t p = 0; p < splitted.size(); ++p) {
        const std::string &entry = splitted[p];
        size_t eq = entry.find('=');
        std::string key = (eq == std::string::npos) ? entry : entry.substr(0, eq);
        if (key == "AF") return (eq == std::string::npos) ? "" : entry.substr(eq + 1);
    }
    return "";
}

// ---------------------------------------------------------------------------
// add_to_dict_snps(line)
// ---------------------------------------------------------------------------

static void add_to_dict_snps(const std::vector<std::string> &line) {
    if (line[3].size() == 1 && line[4].size() == 1) {
        // Case A: biallelic SNP
        std::vector<std::string> list_samples;
        for (size_t p = 9; p < line.size(); ++p) {
            const std::string &i = line[p];
            std::string gt = splitStr(i, ':')[0];
            if (gt.find('1') != std::string::npos) {
                list_samples.push_back(VCFheader[p] + ":" + gt);
            }
        }
        std::string chr_pos_string = currentChr + "," + line[1];
        std::string rsID = line[2];
        std::string af = afValue(line);
        std::string val;
        if (!list_samples.empty()) {
            std::vector<std::string> sorted_s = list_samples;
            std::sort(sorted_s.begin(), sorted_s.end());
            std::string joined;
            for (size_t k = 0; k < sorted_s.size(); ++k) {
                if (k) joined += ",";
                joined += sorted_s[k];
            }
            val = joined + ";" + line[3] + "," + line[4] + ";" + rsID + ";" + af;
        } else {
            val = std::string(";") + line[3] + "," + line[4] + ";" + rsID + ";" + af;
        }
        dictPut(chr_pos_string, std::vector<std::string>(1, val), std::vector<int>(1, 1));
    } else if (line[3].size() == 1) {
        // Case B: multiallelic, ref length 1
        std::vector<std::string> variants = splitStr(line[4], ',');
        std::vector<std::string> snps;
        std::vector<int> values_for_allele_info;
        for (size_t v = 0; v < variants.size(); ++v) {
            if (variants[v].size() == 1) {
                snps.push_back(variants[v]);
                values_for_allele_info.push_back((int)v + 1);
            }
        }
        // dict_of_lists_samples keyed by snp string; python uses a dict so
        // duplicate snp strings collapse. Preserve that: map snp->list.
        std::map<std::string, std::vector<std::string>> dols;
        for (size_t k = 0; k < snps.size(); ++k) dols[snps[k]];  // default-init
        if (!snps.empty()) {
            for (size_t p = 9; p < line.size(); ++p) {
                const std::string &sample = line[p];
                std::string gt = splitStr(sample, ':')[0];
                std::set<int> alleles = gtAltAlleles(gt);
                // No break: a sample carrying two different ALTs (e.g. "2|1")
                // belongs to BOTH of them. Each entry holds its carriers' genotypes in
                // split form (gtForAllele: 1 the entry's own ALT), as a split record's
                // entry and a merged entry do -- enricher.py's add_to_dict_snps.
                for (size_t idx = 0; idx < values_for_allele_info.size(); ++idx) {
                    if (alleles.count(values_for_allele_info[idx])) {
                        dols[snps[idx]].push_back(VCFheader[p] + ":" +
                                                  gtForAllele(gt, values_for_allele_info[idx]));
                    }
                }
            }
            std::string chr_pos_string = currentChr + "," + line[1];
            std::string rsID0 = splitStr(line[2], ',')[0];
            std::vector<std::string> af = splitStr(afValue(line), ',');

            std::vector<std::string> final_entry;
            for (size_t idx = 0; idx < snps.size(); ++idx) {
                std::vector<std::string> &samps = dols[snps[idx]];
                size_t afIdx = (size_t)(values_for_allele_info[idx] - 1);
                // multiallelic site may have fewer AF values than ALT alleles
                std::string afv = (afIdx < af.size()) ? af[afIdx] : "";
                std::string entry;
                if (!samps.empty()) {
                    std::vector<std::string> sorted_s = samps;
                    std::sort(sorted_s.begin(), sorted_s.end());
                    std::string joined;
                    for (size_t k = 0; k < sorted_s.size(); ++k) {
                        if (k) joined += ",";
                        joined += sorted_s[k];
                    }
                    entry = joined + ";" + line[3] + "," + snps[idx] + ";" + rsID0 +
                            ";" + afv;
                } else {
                    entry = std::string(";") + line[3] + "," + snps[idx] + ";" + rsID0 +
                            ";" + afv;
                }
                final_entry.push_back(entry);
            }
            // every entry in split form: its ALT is allele 1 for dictPut's merge
            dictPut(chr_pos_string, final_entry, std::vector<int>(final_entry.size(), 1));
        }
    }
    // len(line[3])>1 : stores nothing
}

// ---------------------------------------------------------------------------
// The fake-indel window (enricher.py's INDEL_LEFT / INDEL_RIGHT): each indel allele's
// contig is INDEL_LEFT reference bases, then the record's REF (from POS on) with the
// ALT written in its place, then INDEL_RIGHT reference bases (fewer where the
// chromosome ends); one N stands between two contigs. A target of L bases (guide +
// PAM) with b DNA bulges holds a record's change only where L + b <= INDEL_LEFT + 2 and
// L + b <= INDEL_RIGHT + 1 (INDEL_LEFT + 1 for a record whose REF and ALT differ at
// their first base), so LEFT = RIGHT = Lmax - 1 + bmax: Lmax = 29 (a 21-nt guide and a
// PAM of 8 nt, the longest the general K=8 index serves), bmax = 2 (the most DNA bulges
// CRISPRme's web offers).
// Up to 2.9.1 the window was 25 / 27 (L + b <= 27). The log's CHR field records each
// window's real start, from which a reader derives the left flank (POS - 1 - start).
// ---------------------------------------------------------------------------

static const long INDEL_LEFT = 30;
static const long INDEL_RIGHT = 30;

// ---------------------------------------------------------------------------
// indel_to_fasta(line, id_indel, start_fake_pos)
// returns updated (id_indel, start_fake_pos)
// ---------------------------------------------------------------------------

static void indel_to_fasta(const std::vector<std::string> &line, long &id_indel,
                           long &start_fake_pos) {
    bool cond = (line[3].size() > 1 || line[4].size() > 1) &&
                line[3].find('<') == std::string::npos;
    if (!cond) return;

    std::vector<std::string> indels;
    std::vector<int> values_for_allele_info;
    if (line[4].find(',') != std::string::npos) {
        std::vector<std::string> splitted = splitStr(line[4], ',');
        for (size_t v = 0; v < splitted.size(); ++v) {
            const std::string &s = splitted[v];
            bool isIndel = (s.size() != line[3].size() ||
                            (s.size() > 1 && line[3].size() > 1)) &&
                           s.find('<') == std::string::npos;
            if (isIndel) {
                indels.push_back(s);
                values_for_allele_info.push_back((int)v + 1);
            }
        }
    } else if (((line[4].size() != line[3].size()) ||
                (line[4].size() > 1 && line[3].size() > 1)) &&
               line[4].find('<') == std::string::npos) {
        indels.push_back(line[4]);
        values_for_allele_info.push_back(1);
    }

    if (indels.empty()) return;

    std::map<std::string, std::vector<std::string>> dols;
    for (size_t k = 0; k < indels.size(); ++k) dols[indels[k]];
    for (size_t p = 9; p < line.size(); ++p) {
        std::string gt = splitStr(line[p], ':')[0];
        std::set<int> alleles = gtAltAlleles(gt);
        // No break: a sample carrying two different ALTs belongs to BOTH.
        for (size_t idx = 0; idx < values_for_allele_info.size(); ++idx) {
            if (alleles.count(values_for_allele_info[idx])) {
                dols[indels[idx]].push_back(VCFheader[p]);
            }
        }
    }

    std::string rsID0 = splitStr(line[2], ',')[0];
    std::vector<std::string> af = splitStr(afValue(line), ',');

    long posv = strtol(line[1].c_str(), NULL, 10);
    long start_position = posv - 1 - INDEL_LEFT;
    long end_position = posv - 1 + (long)line[3].size() + INDEL_RIGHT;
    if (start_position < 0) {
        // Matches enricher.py: the window is anchored with the REF base at index
        // INDEL_LEFT, which a negative start cannot satisfy. substr() would convert the
        // negative offset to size_t and throw, where Python silently yields an empty
        // window -- so the two ports only agree if both skip the site.
        // printf, not cerr: a write to stderr is fatal to the CRISPRme stage that
        // calls this. Matches the file's existing <cstdio> idiom.
        printf("WARNING: skipping indel at %s:%s -- fewer than %ld bases from the chromosome start\n",
               currentChr.c_str(), line[1].c_str(), INDEL_LEFT + 1);
        return;
    }
    std::string sub_fasta = genomeStr.substr(start_position, end_position - start_position);

    for (size_t idx = 0; idx < indels.size(); ++idx) {
        std::vector<std::string> &samps = dols[indels[idx]];
        if (!samps.empty()) {
            std::string indel_info =
                currentChr + "_" + line[1] + "_" + line[3] + "_" + indels[idx];
            // allele_fasta = sub_fasta[0:INDEL_LEFT]
            //               + re.sub(ref, indel, sub_fasta[INDEL_LEFT:], 1, IGNORECASE)
            // Derived from the REFERENCE window each time -- see enricher.py for why
            // assigning back into sub_fasta compounds alleles at a multiallelic site.
            std::string head = sub_fasta.substr(0, INDEL_LEFT);
            std::string tail = sub_fasta.substr(INDEL_LEFT);
            std::string allele_fasta = head + reSubFirstCI(line[3], indels[idx], tail);

            list_fasta_indels.push_back(allele_fasta + "\n" + "N\n");

            std::string refseq = genomeStr.substr(start_position, allele_fasta.size());
            long end_fake_pos = start_fake_pos + (long)allele_fasta.size();

            std::string samplesJoined;
            for (size_t k = 0; k < samps.size(); ++k) {
                if (k) samplesJoined += ",";
                samplesJoined += samps[k];
            }
            std::vector<std::string> row(7);
            row[0] = currentChr + "_" + std::to_string(start_position) + "-" +
                     std::to_string(end_position) + "_" + std::to_string(id_indel);
            row[1] = samplesJoined;
            row[2] = rsID0;
            {
                // multiallelic site may have fewer AF values than ALT alleles
                size_t afIdx = (size_t)(values_for_allele_info[idx] - 1);
                row[3] = (afIdx < af.size()) ? af[afIdx] : "";
            }
            row[4] = indel_info;
            row[5] = std::to_string(start_fake_pos) + "," + std::to_string(end_fake_pos);
            row[6] = refseq;
            log_indels.push_back(row);

            id_indel += 1;
            start_fake_pos = end_fake_pos + 1;
        }
    }
}

// ---------------------------------------------------------------------------
// JSON serialization matching python json.dump defaults:
//   {"k": "v", "k2": "v2"}  -- separators ", " and ": ", insertion order.
// Our values are ASCII; still escape " and \ and control chars per json spec.
// ---------------------------------------------------------------------------

static std::string jsonEscape(const std::string &s) {
    std::string out;
    out.reserve(s.size() + 2);
    for (size_t i = 0; i < s.size(); ++i) {
        unsigned char c = (unsigned char)s[i];
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    char buf[8];
                    snprintf(buf, sizeof(buf), "\\u%04x", c);
                    out += buf;
                } else {
                    out += (char)c;
                }
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// pandas to_csv(index=False, sep='\t') quoting.
// pandas default quoting = QUOTE_MINIMAL: quotes a field only if it contains
// the sep, a quotechar ("), or a line terminator (\n / \r). Commas do NOT
// trigger quoting when sep is a tab. Verified empirically against pandas.
// ---------------------------------------------------------------------------

static std::string csvField(const std::string &s) {
    bool needQuote = false;
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (c == '\t' || c == '"' || c == '\n' || c == '\r') { needQuote = true; break; }
    }
    if (!needQuote) return s;
    std::string out = "\"";
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '"') out += "\"\"";
        else out += s[i];
    }
    out += "\"";
    return out;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char **argv) {
    time_t start_time = time(NULL);
    printf("READING VCF FILE AND CHROMOSOME\n");

    if (argc < 6) {
        fprintf(stderr,
                "Usage: %s <vcf.gz> <chrom.fa> <argv3> <yes|no> <argv5path>\n",
                argv[0]);
        return 1;
    }

    std::string altFile = argv[1];
    std::string genomeFile = argv[2];
    std::string dir_enr_name = std::string(argv[3]) + "_enriched";
    std::string argv4 = argv[4];
    doYes = (argv4 == "yes");

    // vcfName = os.path.basename(argv[5])
    std::string arg5 = argv[5];
    std::string vcfName;
    {
        size_t slash = arg5.find_last_of('/');
        vcfName = (slash == std::string::npos) ? arg5 : arg5.substr(slash + 1);
    }

    // --- read FASTA ---
    FILE *gf = fopen(genomeFile.c_str(), "r");
    if (!gf) {
        fprintf(stderr, "ERROR: cannot open %s\n", genomeFile.c_str());
        return 1;
    }
    // read header line (keep raw, including trailing '\n' like python readline)
    {
        std::string headerLine;
        int c;
        while ((c = fgetc(gf)) != EOF) {
            headerLine.push_back((char)c);
            if (c == '\n') break;
        }
        genomeHeader = headerLine;  // includes '>' and '\n' (if present)
        // currentChr = header[1:].strip()
        std::string h = genomeHeader.substr(genomeHeader.empty() ? 0 : 1);
        currentChr = pyStrip(h);
        if (currentChr.find("chr") == std::string::npos) {
            currentChr = "chr" + currentChr;
        }
    }
    // read rest, remove '\n', uppercase
    {
        std::string rest;
        char buf[1 << 16];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), gf)) > 0) {
            rest.append(buf, n);
        }
        fclose(gf);
        genomeStr.reserve(rest.size());
        for (size_t i = 0; i < rest.size(); ++i) {
            char c = rest[i];
            if (c == '\n') continue;
            if (c >= 'a' && c <= 'z') c -= 32;
            genomeStr.push_back(c);
        }
    }
    // Q is reserved for four-allele SNP positions (see enricher.py)
    {
        size_t q = genomeStr.find('Q');
        if (q != std::string::npos) {
            fprintf(stderr,
                    "ERROR: %s holds 'Q' at %ld: a reference FASTA must not hold Q, which is "
                    "reserved for four-allele SNP positions\n",
                    currentChr.c_str(), (long)q + 1);
            return 1;
        }
    }
    genomeList = genomeStr;  // mutable copy
    snpMask.assign(genomeStr.size(), 0);

    // --- open VCF ---
    GzReader vcf(altFile.c_str());
    // python: VCFheader = inAltFile.readline().split('\t')  (first line, unused;
    // it gets overwritten by the #CHROM search below). Consume first line.
    {
        std::string firstLine;
        vcf.getRawLine(firstLine);  // discarded (matches python overwrite)
    }

    // set up yes-mode structures / directories BEFORE the loop (python does
    // this in original cwd).
    long id_indel = 1;
    long start_fake_pos = 0;
    std::string fakeDir;
    if (doYes) {
        fakeDir = "fake_" + vcfName + "_" + currentChr;
        // race-safe: concurrent workers may create this dir; treat EEXIST as OK.
        if (mkdir(fakeDir.c_str(), 0777) != 0 && errno != EEXIST) {
            fprintf(stderr, "ERROR: cannot create %s\n", fakeDir.c_str());
            return 1;
        }
    }

    // find #CHROM header line
    printf("START ENRICHMENT WITH SNVs AND SVs\n");
    {
        std::string line;
        while (vcf.getRawLine(line)) {
            if (line.find("#CHROM") != std::string::npos) {
                VCFheader = splitStr(pyStrip(line), '\t');
                break;
            }
        }
    }

    // process data lines
    std::string raw;
    while (vcf.getRawLine(raw)) {
        std::vector<std::string> line = splitStr(pyStrip(raw), '\t');
        // issue #143: the per-position variant-info string is ';'-delimited
        // (samples;ref,alt;rsID;AF). The VCF ID column (line[2]) can contain ';'
        // for dbSNP multi-rsID records ("rs1;rs2"), colliding with that delimiter
        // and shifting AF (and every field after) by one -> downstream float(rsID)
        // crash. Normalize the ID field's ';' to ',' so rsID stays one field.
        if (line.size() > 2)
            std::replace(line[2].begin(), line[2].end(), ';', ',');
        if (line[6] != "PASS" && line[6] != ".") continue;  // '.' = no filter applied (e.g. HPRC vcfwave VCFs)
        if (doYes) {
            add_to_dict_snps(line);
            indel_to_fasta(line, id_indel, start_fake_pos);
        }
        SNPsProcess(line);
    }
    if (fourBaseCount)
        printf("NOTE: %ld position(s) of %s are coded Q, the four-allele symbol, which the search "
               "reads as any base (N stays a gap)\n", fourBaseCount, currentChr.c_str());

    // --- chromosomeSave(): os.chdir("./SNPs_genome/") then write ---
    if (chdir("./SNPs_genome") != 0) {
        fprintf(stderr, "ERROR: cannot chdir into ./SNPs_genome\n");
        return 1;
    }
    {
        // race-safe: this shared <dirGenome>_enriched dir is created by every
        // concurrent worker; treat EEXIST as success.
        if (mkdir(dir_enr_name.c_str(), 0777) != 0 && errno != EEXIST) {
            fprintf(stderr, "ERROR: cannot create %s\n", dir_enr_name.c_str());
            return 1;
        }
        // genomeHeader[1:(len-1)] : strip leading '>' and trailing '\n'
        std::string hcore;
        if (genomeHeader.size() >= 2) {
            hcore = genomeHeader.substr(1, genomeHeader.size() - 2);
        }
        std::string outPath = dir_enr_name + "/" + hcore + ".enriched.fa";
        FILE *out = fopen(outPath.c_str(), "w");
        if (!out) {
            fprintf(stderr, "ERROR: cannot write %s\n", outPath.c_str());
            return 1;
        }
        fwrite(genomeHeader.data(), 1, genomeHeader.size(), out);
        fwrite(genomeList.data(), 1, genomeList.size(), out);
        fputc('\n', out);
        fclose(out);
    }

    if (doYes) {
        // dictSave(): writes my_dict_<chr>.json in current cwd (SNPs_genome)
        {
            std::string jpath = "my_dict_" + currentChr + ".json";
            FILE *jf = fopen(jpath.c_str(), "w");
            if (!jf) {
                fprintf(stderr, "ERROR: cannot write %s\n", jpath.c_str());
                return 1;
            }
            std::string js = "{";
            for (size_t k = 0; k < chr_dict_snps.size(); ++k) {
                if (k) js += ", ";
                js += "\"" + jsonEscape(chr_dict_snps[k].first) + "\": \"" +
                      jsonEscape(chr_dict_snps[k].second) + "\"";
            }
            js += "}";
            fwrite(js.data(), 1, js.size(), jf);
            fclose(jf);
        }

        // logIndelsSave():
        //   - fasta_out opened at ORIGINAL cwd path fakeDir/fake<chr>.fa
        //     but we are now in SNPs_genome, so use ../ + fakeDir.
        //   - log<chr>.txt written in current cwd (SNPs_genome).
        {
            std::string fastaPath = "../" + fakeDir + "/fake" + currentChr + ".fa";
            FILE *ff = fopen(fastaPath.c_str(), "w");
            if (!ff) {
                fprintf(stderr, "ERROR: cannot write %s\n", fastaPath.c_str());
                return 1;
            }
            std::string content = ">fake" + currentChr + "\n";
            for (size_t k = 0; k < list_fasta_indels.size(); ++k) {
                content += list_fasta_indels[k];
            }
            fwrite(content.data(), 1, content.size(), ff);
            fclose(ff);
        }
        {
            std::string logPath = "log" + currentChr + ".txt";
            FILE *lf = fopen(logPath.c_str(), "w");
            if (!lf) {
                fprintf(stderr, "ERROR: cannot write %s\n", logPath.c_str());
                return 1;
            }
            const char *cols[7] = {"CHR", "SAMPLES", "rsID", "AF",
                                   "indel", "FAKEPOS", "refseq"};
            std::string out;
            for (int c = 0; c < 7; ++c) {
                if (c) out += "\t";
                out += cols[c];
            }
            out += "\n";
            for (size_t r = 0; r < log_indels.size(); ++r) {
                for (int c = 0; c < 7; ++c) {
                    if (c) out += "\t";
                    out += csvField(log_indels[r][c]);
                }
                out += "\n";
            }
            fwrite(out.data(), 1, out.size(), lf);
            fclose(lf);
        }
    }

    printf("DONE IN  %d  seconds\n", (int)(time(NULL) - start_time));
    return 0;
}

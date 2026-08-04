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

// Maps a mask of ACGT bits to the IUPAC letter reproduced by enricher.py.
static char iupacFromMask(int mask) {
    switch (mask) {
        case 1:  return 'A';
        case 2:  return 'C';
        case 4:  return 'G';
        case 8:  return 'T';
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
        case (1 | 2 | 4 | 8): return 'N';  // A,C,G,T
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

// SNP dict: preserve insertion order (python dict) -> use vector of pairs,
// with a lookup map for "already present" overwrite semantics.
static std::vector<std::pair<std::string, std::string>> chr_dict_snps;
static std::map<std::string, size_t> chr_dict_index;

static void dictSet(const std::string &key, const std::string &val) {
    std::map<std::string, size_t>::iterator it = chr_dict_index.find(key);
    if (it != chr_dict_index.end()) {
        chr_dict_snps[it->second].second = val;  // overwrite, keep position
    } else {
        chr_dict_index[key] = chr_dict_snps.size();
        chr_dict_snps.push_back(std::make_pair(key, val));
    }
}

// fake indel fasta pieces and log rows
static std::vector<std::string> list_fasta_indels;
static std::vector<std::vector<std::string>> log_indels;  // 7 cols each

static bool doYes = false;

// ---------------------------------------------------------------------------
// SNPsProcess(line): mutates genomeList and line[1] -> 0-based.
// ---------------------------------------------------------------------------

static void SNPsProcess(std::vector<std::string> &line) {
    if (line[3].size() > 1) return;  // not a SNP
    long pos0 = strtol(line[1].c_str(), NULL, 10) - 1;
    line[1] = std::to_string(pos0);  // python sets line[1] to 1-based-1

    // build the base SET: reference base + each alt allele of length 1
    int mask = baseBit(genomeStr[pos0]);
    std::vector<std::string> altAlleles = splitStr(pyStrip(line[4]), ',');
    for (size_t k = 0; k < altAlleles.size(); ++k) {
        if (altAlleles[k].size() > 1) continue;
        if (!altAlleles[k].empty()) mask |= baseBit(altAlleles[k][0]);
    }
    char letter = iupacFromMask(mask);
    if (letter != 0) {
        genomeList[pos0] = letter;
    }
    // if letter==0, python loop found no matching value -> genome unchanged.
}

// ---------------------------------------------------------------------------
// add_to_dict_snps(line, pos_AF)
// ---------------------------------------------------------------------------

static void add_to_dict_snps(const std::vector<std::string> &line, int pos_AF) {
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
        std::string af = splitStr(line[7], ';')[pos_AF].substr(3);
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
        dictSet(chr_pos_string, val);
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
                for (size_t idx = 0; idx < values_for_allele_info.size(); ++idx) {
                    std::string vstr = std::to_string(values_for_allele_info[idx]);
                    if (gt.find(vstr) != std::string::npos) {
                        dols[snps[idx]].push_back(VCFheader[p] + ":" + gt);
                        break;
                    }
                }
            }
            std::string chr_pos_string = currentChr + "," + line[1];
            std::string rsID0 = splitStr(line[2], ',')[0];
            std::vector<std::string> af = splitStr(splitStr(line[7], ';')[pos_AF].substr(3), ',');

            std::vector<std::string> final_entry;
            for (size_t idx = 0; idx < snps.size(); ++idx) {
                std::vector<std::string> &samps = dols[snps[idx]];
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
                            ";" + af[values_for_allele_info[idx] - 1];
                } else {
                    entry = std::string(";") + line[3] + "," + snps[idx] + ";" + rsID0 +
                            ";" + af[values_for_allele_info[idx] - 1];
                }
                final_entry.push_back(entry);
            }
            std::string joined;
            for (size_t k = 0; k < final_entry.size(); ++k) {
                if (k) joined += "$";
                joined += final_entry[k];
            }
            dictSet(chr_pos_string, joined);
        }
    }
    // len(line[3])>1 : stores nothing
}

// ---------------------------------------------------------------------------
// indel_to_fasta(line, id_indel, pos_AF, start_fake_pos)
// returns updated (id_indel, start_fake_pos)
// ---------------------------------------------------------------------------

static void indel_to_fasta(const std::vector<std::string> &line, long &id_indel,
                           int pos_AF, long &start_fake_pos) {
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
        for (size_t idx = 0; idx < values_for_allele_info.size(); ++idx) {
            std::string vstr = std::to_string(values_for_allele_info[idx]);
            if (gt.find(vstr) != std::string::npos) {
                dols[indels[idx]].push_back(VCFheader[p]);
                break;
            }
        }
    }

    std::string rsID0 = splitStr(line[2], ',')[0];
    std::vector<std::string> af = splitStr(splitStr(line[7], ';')[pos_AF].substr(3), ',');

    long posv = strtol(line[1].c_str(), NULL, 10);
    long start_position = posv - 26;
    long end_position = posv + 26 + (long)line[3].size();
    std::string sub_fasta = genomeStr.substr(start_position, end_position - start_position);

    for (size_t idx = 0; idx < indels.size(); ++idx) {
        std::vector<std::string> &samps = dols[indels[idx]];
        if (!samps.empty()) {
            std::string indel_info =
                currentChr + "_" + line[1] + "_" + line[3] + "_" + indels[idx];
            // sub_fasta = sub_fasta[0:25] + re.sub(ref, indel, sub_fasta[25:], 1, IGNORECASE)
            std::string head = sub_fasta.substr(0, 25);
            std::string tail = sub_fasta.substr(25);
            sub_fasta = head + reSubFirstCI(line[3], indels[idx], tail);

            list_fasta_indels.push_back(sub_fasta + "\n" + "N\n");

            std::string refseq = genomeStr.substr(start_position, sub_fasta.size());
            long end_fake_pos = start_fake_pos + (long)sub_fasta.size();

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
            row[3] = af[values_for_allele_info[idx] - 1];
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
    genomeList = genomeStr;  // mutable copy

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
    bool first_line = true;
    int pos_AF = 0;
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
        if (first_line) {
            first_line = false;
            std::vector<std::string> splitted = splitStr(line[7], ';');
            for (size_t p = 0; p < splitted.size(); ++p) {
                if (splitted[p].size() >= 2 && splitted[p][0] == 'A' &&
                    splitted[p][1] == 'F') {
                    pos_AF = (int)p;
                    break;
                }
            }
        }
        if (line[6] != "PASS") continue;
        if (doYes) {
            add_to_dict_snps(line, pos_AF);
            indel_to_fasta(line, id_indel, pos_AF, start_fake_pos);
        }
        SNPsProcess(line);
    }

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

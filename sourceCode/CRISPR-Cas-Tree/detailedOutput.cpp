#include <bitset>
#include "convert.h"
#include "detailedOutput.h"
#include <omp.h>

#include <map>
#include <numeric> //for accumulate

inline bool is_true_mismatch(const std::bitset<4>& g, const std::bitset<4>& t) {
    return g.any() && t.any() && (g & t).none();
}

/**This function provides a profiling of the results from the searchOnTST. The results are computed
 * as soon as the target is found on the TST
 * @param guideI: index of the guide 
 * @param guide: the string containing the guide
 * @param target: the string containing the target found
 * @param bulType: int value of the bulge type (0 -> " X ", <0 -> "RNA", >0 -> "DNA")
 * @param mm: number of mismatches allowed
 * @param profiling: an len_guide + mm + 1 x numGuides matrix for profiling.Cells 0-(len_guide-1) (eg 0 - 19) are for counting positional
 * mms, cell 20 for 0MM, cell 21+ for total string count mms
 * devo farmi passare usando & anche le matrici di salvataggio finale
 * @param ext_profiling: matrix for extended profiling. The first dimension selects the guideI, the second dimension selects
 * the matrix to fill w.r.t. the number of mms in the target, the third dimension select the nucleotide (A,C,G,T) and the 
 * fourth dimension selects the position in which there's the mm
 */
void detailedOutputFast(int guideI,
                        std::vector<std::bitset<4>> guide_bit,
                        std::vector<std::bitset<4>> target_bit,
                        std::string &g, std::string &t,
                        int mm, int len_guide,
                        std::vector<std::vector<std::vector<int>>> &profiling,
                        std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling,
                        std::vector<std::string> &vec_guides,
                        std::vector<std::string> &vec_targets,
                        bool pam_at_start)
{
    int thr = omp_get_thread_num();
    int mm_inside_string = 0;
    std::vector<std::pair<std::bitset<4>, int>> save_mm_pos;

    const int L = std::min<int>((int)target_bit.size(), len_guide); // niente PAM

    for (int i = 0; i < L; ++i) {
        if (is_true_mismatch(guide_bit[i], target_bit[i])) {
            ++mm_inside_string;
            profiling[i][guideI][thr]++;
            save_mm_pos.emplace_back(target_bit[i], i);
        }
    }

#ifdef DEBUG
    if (mm_inside_string > mm)
        std::cerr << "Number of mismatches found is > than the one given in input";
#endif

    profiling[len_guide + mm_inside_string][guideI][thr]++;

    for (auto &p : save_mm_pos) {
        if (p.first[0]) ext_profiling[guideI][mm_inside_string][0][p.second][thr]++;
        if (p.first[1]) ext_profiling[guideI][mm_inside_string][1][p.second][thr]++;
        if (p.first[2]) ext_profiling[guideI][mm_inside_string][2][p.second][thr]++;
        if (p.first[3]) ext_profiling[guideI][mm_inside_string][3][p.second][thr]++;
    }

    vec_targets.push_back(t);
    vec_guides.push_back(g);
}

void detailedOutputFastBulgeDNA(int guideI,
                                std::vector<std::bitset<4>> guide_bit,
                                std::vector<std::bitset<4>> target_bit,
                                std::string &g, std::string &t,
                                int mm, int max_bulges, int len_guide, int bD, int bulDNA,
                                std::vector<std::vector<std::vector<int>>> &profiling_dna,
                                std::vector<std::vector<std::vector<int>>> &profiling_dna_mm,
                                std::vector<std::string> &vec_guides,
                                std::vector<std::string> &vec_targets,
                                bool pam_at_start,
                                std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_dna,
                                std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling)
{
    int thr = omp_get_thread_num();
    int mm_inside_string = 0;
    int bul_inside_string = 0; 
    std::vector<int> pos_bulges;
    std::vector<std::pair<std::bitset<4>, int>> save_mm_pos;

    const int L = std::min<int>((int)guide_bit.size(), len_guide + bulDNA - bD);

    for (int i = 0; i < L; ++i) {
        if ((guide_bit[i] & target_bit[i]).any()) continue; 

        if (guide_bit[i].none()) {
            if (i >= (len_guide + bul_inside_string)) {
                profiling_dna[i - bul_inside_string - 1][guideI][thr]++;
                pos_bulges.push_back(i - bul_inside_string - 1);
            } else {
                profiling_dna[i - bul_inside_string][guideI][thr]++;
                pos_bulges.push_back(i - bul_inside_string);
            }
            ++bul_inside_string;
            continue;
        }

        if (is_true_mismatch(guide_bit[i], target_bit[i])) {
            ++mm_inside_string;
            profiling_dna_mm[i - bul_inside_string][guideI][thr]++;
            save_mm_pos.emplace_back(target_bit[i], i - bul_inside_string);
        }
    }

    profiling_dna_mm[len_guide + bulDNA + mm_inside_string * max_bulges + (bul_inside_string - 1)][guideI][thr]++;

    for (int p : pos_bulges)
        ext_profiling_dna[guideI][mm_inside_string][p][thr]++;

    for (auto &p : save_mm_pos) {
        if (p.first[0]) ext_profiling[guideI][mm_inside_string][0][p.second][thr]++;
        if (p.first[1]) ext_profiling[guideI][mm_inside_string][1][p.second][thr]++;
        if (p.first[2]) ext_profiling[guideI][mm_inside_string][2][p.second][thr]++;
        if (p.first[3]) ext_profiling[guideI][mm_inside_string][3][p.second][thr]++;
    }

    vec_guides.push_back(g);
    vec_targets.push_back(t);
}

void detailedOutputFastBulgeRNA(int guideI,
                                std::vector<std::bitset<4>> guide_bit,
                                std::vector<std::bitset<4>> target_bit,
                                std::string &g, std::string &t,
                                int mm, int max_bulges, int len_guide, int bD, int bulDNA,
                                std::vector<std::vector<std::vector<int>>> &profiling_rna,
                                std::vector<std::vector<std::vector<int>>> &profiling_rna_mm,
                                std::vector<std::string> &vec_guides,
                                std::vector<std::string> &vec_targets,
                                bool pam_at_start,
                                std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_rna,
                                std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling)
{
    int thr = omp_get_thread_num();
    int mm_inside_string = 0;
    int bul_inside_string = 0; 
    std::vector<int> pos_bulges;
    std::vector<std::pair<std::bitset<4>, int>> save_mm_pos;

    const int L = std::min<int>((int)guide_bit.size(), len_guide);

    for (int i = 0; i < L; ++i) {
        if ((guide_bit[i] & target_bit[i]).any()) continue; // match

        if (target_bit[i].none()) {
            profiling_rna[i][guideI][thr]++;
            pos_bulges.push_back(i);
            ++bul_inside_string;
            continue;
        }

        // mismatch
        if (is_true_mismatch(guide_bit[i], target_bit[i])) {
            ++mm_inside_string;
            profiling_rna_mm[i][guideI][thr]++;
            save_mm_pos.emplace_back(target_bit[i], i);
        }
    }

    profiling_rna_mm[len_guide + mm_inside_string * max_bulges + (bul_inside_string - 1)][guideI][thr]++;

    for (int p : pos_bulges)
        ext_profiling_rna[guideI][mm_inside_string][p][thr]++;

    for (auto &p : save_mm_pos) {
        if (p.first[0]) ext_profiling[guideI][mm_inside_string][0][p.second][thr]++;
        if (p.first[1]) ext_profiling[guideI][mm_inside_string][1][p.second][thr]++;
        if (p.first[2]) ext_profiling[guideI][mm_inside_string][2][p.second][thr]++;
        if (p.first[3]) ext_profiling[guideI][mm_inside_string][3][p.second][thr]++;
    }

    vec_guides.push_back(g);
    vec_targets.push_back(t);
}


// helper: flat index for the joint (mm, bulgeDNA, bulgeRNA) bin
inline int both_bin_index(int mmv, int bd, int br,
                          int len_guide, int bulDNA, int max_bulges) {
    const int base   = len_guide + bulDNA;
    const int stride = (max_bulges + 1) * (max_bulges + 1);
    return base + mmv * stride + bd * (max_bulges + 1) + br;
}

void detailedOutputFastBulgeBoth(int guideI,
                                 std::vector<std::bitset<4>> guide_bit,
                                 std::vector<std::bitset<4>> target_bit,
                                 std::string &g, std::string &t,
                                 int mm, int max_bulges, int len_guide, int bD, int bulDNA,
                                 std::vector<std::vector<std::vector<int>>> &profiling_dna,
                                 std::vector<std::vector<std::vector<int>>> &profiling_dna_mm,
                                 std::vector<std::vector<std::vector<int>>> &profiling_rna,
                                 std::vector<std::vector<std::vector<int>>> &profiling_rna_mm,
                                 std::vector<std::vector<std::vector<int>>> &profiling_both_mm,
                                 std::vector<std::string> &vec_guides,
                                 std::vector<std::string> &vec_targets,
                                 bool pam_at_start,
                                 std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_dna,
                                 std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_rna,
                                 std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling)
{
    int thr = omp_get_thread_num();
    int mm_inside_string = 0;
    int bul_inside_dna   = 0; // gap in guide
    int bul_inside_rna   = 0; // gap in target

    std::vector<int> pos_bulges_dna, pos_bulges_rna;
    std::vector<std::pair<std::bitset<4>, int>> save_mm_pos;

    const int L = std::min<int>((int)guide_bit.size(), len_guide + bulDNA - bD);

    for (int i = 0; i < L; ++i) {
        if ((guide_bit[i] & target_bit[i]).any()) continue; // match

        // bulge DNA
        if (guide_bit[i].none()) {
            if (i >= (len_guide + bul_inside_dna)) {
                profiling_dna[i - bul_inside_dna - 1][guideI][thr]++;
                pos_bulges_dna.push_back(i - bul_inside_dna - 1);
            } else {
                profiling_dna[i - bul_inside_dna][guideI][thr]++;
                pos_bulges_dna.push_back(i - bul_inside_dna);
            }
            ++bul_inside_dna;
            continue;
        }

        // bulge RNA
        if (i < len_guide && target_bit[i].none()) {
            profiling_rna[i][guideI][thr]++;
            pos_bulges_rna.push_back(i);
            ++bul_inside_rna;
            continue;
        }

        //  mismatch
        if (is_true_mismatch(guide_bit[i], target_bit[i])) {
            ++mm_inside_string;
            profiling_dna_mm[i - bul_inside_dna][guideI][thr]++;
            profiling_rna_mm[i][guideI][thr]++;
            save_mm_pos.emplace_back(target_bit[i], i - bul_inside_dna);
        }
    }

    profiling_dna_mm[len_guide + bulDNA + mm_inside_string * max_bulges + (bul_inside_dna - 1)][guideI][thr]++;
    profiling_rna_mm[len_guide + mm_inside_string * max_bulges + (bul_inside_rna - 1)][guideI][thr]++;

    if (mm_inside_string <= mm && bul_inside_dna <= max_bulges && bul_inside_rna <= max_bulges) {
        const int idx_joint = both_bin_index(mm_inside_string, bul_inside_dna, bul_inside_rna,
                                             len_guide, bulDNA, max_bulges);
        profiling_both_mm[idx_joint][guideI][thr]++;
    }

    for (int p : pos_bulges_dna) ext_profiling_dna[guideI][mm_inside_string][p][thr]++;
    for (int p : pos_bulges_rna) ext_profiling_rna[guideI][mm_inside_string][p][thr]++;

    for (auto &p : save_mm_pos) {
        if (p.first[0]) ext_profiling[guideI][mm_inside_string][0][p.second][thr]++;
        if (p.first[1]) ext_profiling[guideI][mm_inside_string][1][p.second][thr]++;
        if (p.first[2]) ext_profiling[guideI][mm_inside_string][2][p.second][thr]++;
        if (p.first[3]) ext_profiling[guideI][mm_inside_string][3][p.second][thr]++;
    }

    vec_guides.push_back(g);
    vec_targets.push_back(t);
}

/**
 * Extended version: write out profiles for a guide, reducing per-thread planes and
 * using a joint (mm, bulgeDNA, bulgeRNA) matrix to avoid double counting
 * when a target has both DNA and RNA bulges.
 *
 * Parameters are as in your original saveProfileGuide, plus:
 *  - profiling_both_mm: [bins][guide][thr], bins laid out by both_bin_index(...)
 */
void saveProfileGuide(std::string guide, int guideI, int mism, int max_bulges, int len_guide, int bulDNA,
                          std::vector<std::vector<std::vector<int>>> &profiling, // base (pos & mm-only summary)
                          std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling,
                          std::vector<std::vector<std::vector<int>>> &profiling_dna,
                          std::vector<std::vector<std::vector<int>>> &profiling_dna_mm,
                          std::vector<std::vector<std::vector<int>>> &profiling_rna,
                          std::vector<std::vector<std::vector<int>>> &profiling_rna_mm,
                          std::vector<std::vector<std::vector<int>>> &profiling_both_mm, // NEW joint bins
                          std::ofstream &fileprofiling, std::ofstream &file_ext_profiling,
                          std::ofstream &file_profiling_dna, std::ofstream &file_profiling_rna,
                          std::ofstream &file_profiling_complete,
                          int num_thr, int pam_size, bool pam_at_start,
                          std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_dna,
                          std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_rna)
{
    // 1) Reduce per-thread planes
    // base mismatch pos & mm-summary
    for (int i = 0; i < (int)profiling.size(); ++i)
        for (int th = 1; th < num_thr; ++th)
            profiling[i][guideI][0] += profiling[i][guideI][th];

    // RNA pos + summary
    for (int i = 0; i < (int)profiling_rna.size(); ++i)
        for (int th = 1; th < num_thr; ++th) {
            profiling_rna[i][guideI][0]    += profiling_rna[i][guideI][th];
            profiling_rna_mm[i][guideI][0] += profiling_rna_mm[i][guideI][th];
        }

    // DNA pos + summary
    for (int i = 0; i < (int)profiling_dna_mm.size(); ++i)
        for (int th = 1; th < num_thr; ++th) {
            profiling_dna[i][guideI][0]    += profiling_dna[i][guideI][th];
            profiling_dna_mm[i][guideI][0] += profiling_dna_mm[i][guideI][th];
        }

    // NEW: joint (mm,bd,br) summary
    for (int i = 0; i < (int)profiling_both_mm.size(); ++i)
        for (int th = 1; th < num_thr; ++th)
            profiling_both_mm[i][guideI][0] += profiling_both_mm[i][guideI][th];

    for (int mmv = 0; mmv <= mism; ++mmv)
        for (int nuc = 0; nuc < 4; ++nuc)
            for (int pos = 0; pos < len_guide; ++pos)
                for (int th = 1; th < num_thr; ++th)
                    ext_profiling[guideI][mmv][nuc][pos][0] += ext_profiling[guideI][mmv][nuc][pos][th];

    for (int mmv = 0; mmv <= mism; ++mmv)
        for (int pos = 0; pos < len_guide; ++pos)
            for (int th = 1; th < num_thr; ++th) {
                ext_profiling_dna[guideI][mmv][pos][0] += ext_profiling_dna[guideI][mmv][pos][th];
                ext_profiling_rna[guideI][mmv][pos][0] += ext_profiling_rna[guideI][mmv][pos][th];
            }

    // 2) Label (guide + PAM Ns only for printing)
    if (pam_at_start) guide.insert(0, std::string(pam_size, 'N'));
    else              guide.append(std::string(pam_size, 'N'));

    std::vector<int> bp_sum_complete(len_guide + bulDNA + (mism+1)*(max_bulges+1)*(max_bulges+1), 0);
    int ont_complete = 0, offt_complete = 0;

    fileprofiling << guide << "\t";
    // per-position mismatch counts
    int sum_bp = 0;
    for (int i = 0; i < len_guide; ++i) {
        fileprofiling << profiling[i][guideI][0] << "\t";
        bp_sum_complete[i] += profiling[i][guideI][0];
        sum_bp += profiling[i][guideI][0];
    }
    fileprofiling << "\t";

    fileprofiling << profiling[len_guide][guideI][0] << "\t";
    ont_complete += profiling[len_guide][guideI][0];

    int sum_offt_base = 0;
    for (int i = len_guide + 1; i < (int)profiling.size(); ++i)
        sum_offt_base += profiling[i][guideI][0];
    offt_complete += sum_offt_base;

    fileprofiling << sum_offt_base << "\t"
                  << (sum_offt_base ? (sum_bp / (double)sum_offt_base) : 0.0) << "\t\t";

    for (int i = len_guide; i < (int)profiling.size(); ++i) {
        fileprofiling << profiling[i][guideI][0] << "\t";
        bp_sum_complete[i + bulDNA] += profiling[i][guideI][0];
    }
    fileprofiling << "\n";

    // 4) DNA profiling
    file_profiling_dna << guide << "\t";
    sum_bp = 0;
    for (int i = 0; i < len_guide; ++i) {
        file_profiling_dna << profiling_dna_mm[i][guideI][0] << "("
                           << profiling_dna[i][guideI][0] << ")\t";
        bp_sum_complete[i] += profiling_dna_mm[i][guideI][0] + profiling_dna[i][guideI][0];
        sum_bp += profiling_dna_mm[i][guideI][0];
    }
    file_profiling_dna << "\t";

    int total_ont_dna = 0; // sum of 0MM with j bulges (DNA channel)
    for (int j = 0; j < max_bulges; ++j)
        total_ont_dna += profiling_dna_mm[len_guide + bulDNA + j][guideI][0];
    ont_complete += total_ont_dna;
    file_profiling_dna << total_ont_dna << "\t";

    int sum_offt_dna = 0;
    for (int i = len_guide + bulDNA + max_bulges; i < (int)profiling_dna_mm.size(); ++i)
        sum_offt_dna += profiling_dna_mm[i][guideI][0];
    offt_complete += sum_offt_dna;

    file_profiling_dna << sum_offt_dna << "\t"
                       << (sum_offt_dna ? (sum_bp / (double)sum_offt_dna) : 0.0) << "\t\t";

    for (int i = len_guide + bulDNA; i < (int)profiling_dna_mm.size(); ++i)
        file_profiling_dna << profiling_dna_mm[i][guideI][0] << "\t";
    file_profiling_dna << "\n";

    // 5) RNA profiling
    file_profiling_rna << guide << "\t";
    sum_bp = 0;
    for (int i = 0; i < len_guide; ++i) {
        file_profiling_rna << profiling_rna_mm[i][guideI][0] << "("
                           << profiling_rna[i][guideI][0] << ")\t";
        bp_sum_complete[i] += profiling_rna_mm[i][guideI][0] + profiling_rna[i][guideI][0];
        sum_bp += profiling_rna_mm[i][guideI][0];
    }
    file_profiling_rna << "\t";

    int total_ont_rna = 0; // sum of 0MM with j bulges (RNA channel)
    for (int j = 0; j < max_bulges; ++j)
        total_ont_rna += profiling_rna_mm[len_guide + j][guideI][0];
    ont_complete += total_ont_rna;
    file_profiling_rna << total_ont_rna << "\t";

    int sum_offt_rna = 0;
    for (int i = len_guide + max_bulges; i < (int)profiling_rna_mm.size(); ++i)
        sum_offt_rna += profiling_rna_mm[i][guideI][0];
    offt_complete += sum_offt_rna;

    file_profiling_rna << sum_offt_rna << "\t"
                       << (sum_offt_rna ? (sum_bp / (double)sum_offt_rna) : 0.0) << "\t\t";

    for (int i = len_guide; i < (int)profiling_rna_mm.size(); ++i)
        file_profiling_rna << profiling_rna_mm[i][guideI][0] << "\t";
    file_profiling_rna << "\n";

    // 6) COMPLETE profiling (now using the joint matrix to avoid double counting)
    file_profiling_complete << guide << "\t";

    for (int i = 0; i < len_guide; ++i)
        file_profiling_complete << bp_sum_complete[i] << "\t";

    file_profiling_complete << "\t";

    int pos_bp_sum = len_guide + bulDNA;
    for (int mmv = 0; mmv <= mism; ++mmv) {
        for (int bd = 0; bd <= max_bulges; ++bd) {
            for (int br = 0; br <= max_bulges; ++br) {
                int idx = both_bin_index(mmv, bd, br, len_guide, bulDNA, max_bulges);
                bp_sum_complete[ idx ] += profiling_both_mm[ idx ][guideI][0];
            }
        }
    }

    int sum_bp_complete = 0;
    for (int i = 0; i < len_guide; ++i) sum_bp_complete += bp_sum_complete[i];

    file_profiling_complete << ont_complete << "\t"
                            << offt_complete << "\t"
                            << (offt_complete ? (sum_bp_complete / (double)offt_complete) : 0.0) << "\t\t";

    for (int i = pos_bp_sum; i < (int)bp_sum_complete.size(); ++i)
        file_profiling_complete << bp_sum_complete[i] << "\t";
    file_profiling_complete << "\n";

    // 7) EXT profiling 
    file_ext_profiling << ">" << guide << "\n\t";
    for (int j = 0; j < len_guide; ++j) file_ext_profiling << "BP\t";
    file_ext_profiling << "TARGETS\n";

    int scan_idx = pos_bp_sum;
    for (int mmv = 0; mmv <= mism; ++mmv)
    {
        file_ext_profiling << mmv << " MISMATCHES\t";
        for (int j = 0; j < len_guide; ++j) {
            int sum = 0;
            for (int k = 0; k < 4; ++k) sum += ext_profiling[guideI][mmv][k][j][0];
            sum += ext_profiling_dna[guideI][mmv][j][0] + ext_profiling_rna[guideI][mmv][j][0];
            file_ext_profiling << sum << "\t";
        }
        // total targets for this mm from the complete (joint) block:
        file_ext_profiling << bp_sum_complete[scan_idx] << "\n";
        ++scan_idx;

        file_ext_profiling << "NUCLEOTIDE A\t";
        for (int j = 0; j < len_guide; ++j) file_ext_profiling << ext_profiling[guideI][mmv][0][j][0] << "\t";
        file_ext_profiling << "\n";

        file_ext_profiling << "NUCLEOTIDE C\t";
        for (int j = 0; j < len_guide; ++j) file_ext_profiling << ext_profiling[guideI][mmv][1][j][0] << "\t";
        file_ext_profiling << "\n";

        file_ext_profiling << "NUCLEOTIDE G\t";
        for (int j = 0; j < len_guide; ++j) file_ext_profiling << ext_profiling[guideI][mmv][2][j][0] << "\t";
        file_ext_profiling << "\n";

        file_ext_profiling << "NUCLEOTIDE T\t";
        for (int j = 0; j < len_guide; ++j) file_ext_profiling << ext_profiling[guideI][mmv][3][j][0] << "\t";
        file_ext_profiling << "\n";

        file_ext_profiling << "Bulge DNA\t";
        for (int j = 0; j < len_guide; ++j) file_ext_profiling << ext_profiling_dna[guideI][mmv][j][0] << "\t";
        file_ext_profiling << "\n";

        file_ext_profiling << "Bulge RNA\t";
        for (int j = 0; j < len_guide; ++j) file_ext_profiling << ext_profiling_rna[guideI][mmv][j][0] << "\t";
        file_ext_profiling << "\n\n";
    }
}
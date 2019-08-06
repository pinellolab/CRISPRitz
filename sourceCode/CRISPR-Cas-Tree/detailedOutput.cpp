#include <bitset>
#include "convert.h"
#include "detailedOutput.h"
#include <omp.h>

#include <map>
#include <numeric> //for accumulate

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
void detailedOutputFast(int guideI, std::vector<std::bitset<4>> guide_bit, std::vector<std::bitset<4>> target_bit, std::string &g, std::string &t,
                        int bulType, int mm, int len_guide, std::vector<std::vector<std::vector<int>>> &profiling, std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling,
                        std::vector<std::string> &vec_guides, std::vector<std::string> &vec_targets, bool pam_at_start)
{
    // std::cout << "profiling size: " << profiling.size() << std::endl;
    // std::cout << "target: " << t << std::endl;
    int thr = omp_get_thread_num();
    int mm_inside_string = 0;
    std::vector<std::pair<std::bitset<4>, int>> save_mm_pos; //array for the pair (base mm in bit, pos)
    
    for (int i = 0; i < target_bit.size(); i++)
    {
        if (i < len_guide) //check for mismatches only in guide char, not pam
            if ((guide_bit[i] & target_bit[i]) == 0)
            { //mismatch case
                mm_inside_string++;
                profiling[i][guideI][thr]++;
                save_mm_pos.push_back(std::make_pair(target_bit[i], i));
            }
    }
    if (mm_inside_string > mm)
        std::cerr << "Number of mismatches found is > that the one given in input";
    profiling[len_guide + mm_inside_string][guideI][thr]++;
    
    //extended matrix with save_mm_pos
    for (int i = 0; i < save_mm_pos.size(); i++)
    {
        if (save_mm_pos[i].first[0])
        {
            ext_profiling[guideI][mm_inside_string][0][save_mm_pos[i].second][thr]++; //a
        }
        if (save_mm_pos[i].first[1])
        {
            ext_profiling[guideI][mm_inside_string][1][save_mm_pos[i].second][thr]++; //c
        }
        if (save_mm_pos[i].first[2])
        {
            ext_profiling[guideI][mm_inside_string][2][save_mm_pos[i].second][thr]++; //g
        }
        if (save_mm_pos[i].first[3])
        {
            ext_profiling[guideI][mm_inside_string][3][save_mm_pos[i].second][thr]++; //t
        }
    }

    vec_targets.push_back(t);
    vec_guides.push_back(g);
}

void detailedOutputFastBulgeDNA(int guideI, std::vector<std::bitset<4>> guide_bit, std::vector<std::bitset<4>> target_bit, std::string &g, std::string &t,
                                int bulType, int mm, int max_bulges, int len_guide, int bD, int bulDNA, std::vector<std::vector<std::vector<int>>> &profiling_dna,
                                std::vector<std::vector<std::vector<int>>> &profiling_dna_mm,
                                std::vector<std::string> &vec_guides, std::vector<std::string> &vec_targets, bool pam_at_start,
                                std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_dna,
                                std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling)
{

    int thr = omp_get_thread_num();
    int mm_inside_string = 0;
    int bul_inside_string = 0;
    std::vector<int> pos_bulges;
    std::vector<std::pair<std::bitset<4>, int>> save_mm_pos; //array for the pair (base mm in bit, pos)
    for (int i = 0; i < guide_bit.size(); i++)
    {

        if (i < len_guide + bulDNA - bD)
            if ((guide_bit[i] & target_bit[i]) == 0)
            { //mismatch case
                if (guide_bit[i].none())
                {
                    if ( i >= (len_guide + bul_inside_string)){
                        profiling_dna[i - bul_inside_string - 1][guideI][thr]++;
                        pos_bulges.push_back(i - bul_inside_string - 1);
                        bul_inside_string++;
                        continue;
                    }

                    profiling_dna[i - bul_inside_string][guideI][thr]++;
                    pos_bulges.push_back(i - bul_inside_string);
                    bul_inside_string++;
                    continue;
                }
                mm_inside_string++;
                profiling_dna_mm[i - bul_inside_string][guideI][thr]++;
                save_mm_pos.push_back(std::make_pair(target_bit[i], i - bul_inside_string)); // controllo i
            }
    }
    profiling_dna_mm[len_guide + bulDNA + mm_inside_string * max_bulges + (bul_inside_string - 1)][guideI][thr]++;
    for (auto p : pos_bulges)
        ext_profiling_dna[guideI][mm_inside_string][p][thr]++;
    
    //extended matrix with save_mm_pos
    for (int i = 0; i < save_mm_pos.size(); i++)
    {
        if (save_mm_pos[i].first[0])
        {
            ext_profiling[guideI][mm_inside_string][0][save_mm_pos[i].second][thr]++; //a
        }
        if (save_mm_pos[i].first[1])
        {
            ext_profiling[guideI][mm_inside_string][1][save_mm_pos[i].second][thr]++; //c
        }
        if (save_mm_pos[i].first[2])
        {
            ext_profiling[guideI][mm_inside_string][2][save_mm_pos[i].second][thr]++; //g
        }
        if (save_mm_pos[i].first[3])
        {
            ext_profiling[guideI][mm_inside_string][3][save_mm_pos[i].second][thr]++; //t
        }
    }
    vec_guides.push_back(g);
    vec_targets.push_back(t);
}

void detailedOutputFastBulgeRNA(int guideI, std::vector<std::bitset<4>> guide_bit, std::vector<std::bitset<4>> target_bit, std::string &g, std::string &t,
                                int bulType, int mm, int max_bulges, int len_guide, int bD, int bulDNA, std::vector<std::vector<std::vector<int>>> &profiling_rna,
                                std::vector<std::vector<std::vector<int>>> &profiling_rna_mm,
                                std::vector<std::string> &vec_guides, std::vector<std::string> &vec_targets, bool pam_at_start,
                                std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_rna,
                                std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling)
{

    int thr = omp_get_thread_num();
    int mm_inside_string = 0;
    int bul_inside_string = 0;
    std::vector<int> pos_bulges;
    std::vector<std::pair<std::bitset<4>, int>> save_mm_pos; //array for the pair (base mm in bit, pos)

    for (int i = 0; i < guide_bit.size(); i++)
    {
        if (i < len_guide)
            if ((guide_bit[i] & target_bit[i]) == 0)
            { //mismatch case
                if (target_bit[i].none())
                {
                    profiling_rna[i][guideI][thr]++;
                    bul_inside_string++;
                    pos_bulges.push_back(i);
                    continue;
                }
                mm_inside_string++;
                profiling_rna_mm[i][guideI][thr]++;
                save_mm_pos.push_back(std::make_pair(target_bit[i], i));
            }
    }
    profiling_rna_mm[len_guide + mm_inside_string * max_bulges + (bul_inside_string - 1)][guideI][thr]++;
    for (auto p : pos_bulges)
        ext_profiling_rna[guideI][mm_inside_string][p][thr]++;

    //extended matrix with save_mm_pos
    for (int i = 0; i < save_mm_pos.size(); i++)
    {
        if (save_mm_pos[i].first[0])
        {
            ext_profiling[guideI][mm_inside_string][0][save_mm_pos[i].second][thr]++; //a
        }
        if (save_mm_pos[i].first[1])
        {
            ext_profiling[guideI][mm_inside_string][1][save_mm_pos[i].second][thr]++; //c
        }
        if (save_mm_pos[i].first[2])
        {
            ext_profiling[guideI][mm_inside_string][2][save_mm_pos[i].second][thr]++; //g
        }
        if (save_mm_pos[i].first[3])
        {
            ext_profiling[guideI][mm_inside_string][3][save_mm_pos[i].second][thr]++; //t
        }
    }

    vec_guides.push_back(g);
    vec_targets.push_back(t);
}
void saveProfileGuide(std::string guide, int guideI, int mism, int max_bulges, int len_guide, int bulDNA, std::vector<std::vector<std::vector<int>>> &profiling,
                      std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling,
                      std::vector<std::vector<std::vector<int>>> &profiling_dna, std::vector<std::vector<std::vector<int>>> &profiling_dna_mm,
                      std::vector<std::vector<std::vector<int>>> &profiling_rna, std::vector<std::vector<std::vector<int>>> &profiling_rna_mm,
                      std::ofstream &fileprofiling, std::ofstream &file_ext_profiling,
                      std::ofstream &file_profiling_dna, std::ofstream &file_profiling_rna, std::ofstream &file_profiling_complete, int num_thr, int pam_size, bool pam_at_start,
                      std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_dna,
                      std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_rna)
{
    //sum the layer of the matrix
    for (int i = 0; i < profiling.size(); i++)
    {
        for (int j = 1; j < num_thr; j++)
        {
            profiling[i][guideI][0] += profiling[i][guideI][j];
        }
    }

    for (int i = 0; i < profiling_rna.size(); i++){
        for (int j = 1; j < num_thr; j++){
            profiling_rna[i][guideI][0] += profiling_rna[i][guideI][j];
            profiling_rna_mm[i][guideI][0] += profiling_rna_mm[i][guideI][j];
        }
    }

    for (int i = 0; i < profiling_dna_mm.size(); i++)
    {
        for (int j = 1; j < num_thr; j++)
        {
            profiling_dna[i][guideI][0] += profiling_dna[i][guideI][j];
            profiling_dna_mm[i][guideI][0] += profiling_dna_mm[i][guideI][j];
        }
    }
    


    

    for (int mm = 0; mm <= mism; mm++)
    {
        for (int nuc = 0; nuc < 4; nuc++){
            for (int base = 0; base < len_guide; base++){
                for (int j = 1; j < num_thr; j++){
                    ext_profiling[guideI][mm][nuc][base][0] += ext_profiling[guideI][mm][nuc][base][j];
                }
            }
        }
    }

    for (int mm = 0; mm <= mism; mm++){
        for (int pos = 0; pos < len_guide; pos ++){
            for (int j = 1; j < num_thr; j++){
                ext_profiling_dna[guideI][mm][pos][0] += ext_profiling_dna[guideI][mm][pos][j];
                ext_profiling_rna[guideI][mm][pos][0] += ext_profiling_rna[guideI][mm][pos][j];
            }
        }
    }




    //For the complete profile
    std::vector<int> bp_sum_complete(len_guide + bulDNA + mism + 1, 0);
    int ont_complete = 0;
    int offt_complete = 0;
    //save profiling
    if (pam_at_start)
    {
        for (int i = 0; i < pam_size; i++)
            guide.insert(0,"N");
    }
    else
    {
        for (int i = 0; i < pam_size; i++)
            guide+="N";
    }
    fileprofiling << guide << "\t";
    for (int i = 0; i < len_guide; i++)
    {
        fileprofiling << profiling[i][guideI][0] << "\t"; //write BP columns
        bp_sum_complete[i] += profiling[i][guideI][0]; //for complete profiling
    }
    fileprofiling << "\t";
    fileprofiling << profiling[len_guide][guideI][0] << "\t"; //write column of ONT, take data from col with number of 0mms
    ont_complete += profiling[len_guide][guideI][0];
    int sum = 0;
    for (int i = (len_guide + 1); i < profiling.size(); i++)
    {       //sum the column 1MM,2MM...
        sum += profiling[i][guideI][0];
    }
    offt_complete += sum;
    int sum_bp = 0;
    for (int i = 0; i < len_guide; i++){
        sum_bp += profiling[i][guideI][0];
    }
    fileprofiling << sum << "\t" << sum_bp / (double)sum << "\t\t"; //write Total number of mms / total OFFT
    for (int i = len_guide; i < profiling.size(); i++)
    { //write 0MM, 1MM etc
        fileprofiling << profiling[i][guideI][0] << "\t";
        bp_sum_complete[i + bulDNA] += profiling[i][guideI][0];
    }
    fileprofiling << "\n";

    //profiling dna
    file_profiling_dna << guide << "\t";
    for (int i = 0; i < (len_guide); i++) //+bulDNA
    {
        file_profiling_dna << profiling_dna_mm[i][guideI][0] << "(" << profiling_dna[i][guideI][0] << ")"
                           << "\t"; //write BP columns
        bp_sum_complete[i] += profiling_dna_mm[i][guideI][0] + profiling_dna[i][guideI][0];
    }
    file_profiling_dna << "\t";
    int total_ont = 0;
    for (int i = 0; i < max_bulges; i++){
        total_ont += profiling_dna_mm[len_guide + bulDNA + i][guideI][0]; //+ profiling_dna_mm[len_guide + bulDNA + 1][guideI][0]; // sum cols 0mm1 0mm2
    }
    ont_complete += total_ont;
    file_profiling_dna << total_ont << "\t";
    sum = 0;
    for (int i = (len_guide + bulDNA + max_bulges); i < profiling_dna_mm.size(); i++)
    { //+2 (+max_bulges) because i have to skip the columns 0mm1 0mm2
        sum += profiling_dna_mm[i][guideI][0];
    }
    offt_complete += sum;
    sum_bp = 0;
    for (int i = 0; i < len_guide; i++){ //+bulDNA
        sum_bp += profiling_dna_mm[i][guideI][0];
    }
    file_profiling_dna << sum << "\t" << sum_bp / (double)sum << "\t\t";
    for (int i = (len_guide + bulDNA); i < profiling_dna_mm.size(); i++)
    {
        file_profiling_dna << profiling_dna_mm[i][guideI][0] << "\t";
    }
    file_profiling_dna << "\n";

    //profiling rna
    file_profiling_rna << guide << "\t";
    for (int i = 0; i < len_guide; i++)
    {
        file_profiling_rna << profiling_rna_mm[i][guideI][0] << "(" << profiling_rna[i][guideI][0] << ")"
                           << "\t"; //write BP columns: number of mm in that position (number of bulges in that position)
        bp_sum_complete[i] += profiling_rna_mm[i][guideI][0] + profiling_rna[i][guideI][0];
    }
    file_profiling_rna << "\t";
    total_ont = 0;
    for(int i = 0; i < max_bulges; i++){
        total_ont += profiling_rna_mm[len_guide + i][guideI][0];// + profiling_rna_mm[len_guide + 1][guideI][0];
    }

    ont_complete += total_ont;
    file_profiling_rna << total_ont << "\t";
    sum = 0;
    for (int i = (len_guide + max_bulges); i < profiling_rna_mm.size(); i++)
    {
        sum += profiling_rna_mm[i][guideI][0];
    }
    offt_complete += sum;
    sum_bp = 0;
    for (int i = 0; i < len_guide; i++){
        sum_bp += profiling_rna_mm[i][guideI][0];
    }
    file_profiling_rna << sum << "\t" << sum_bp / (double)sum << "\t\t";
    for (int i = (len_guide); i < profiling_rna_mm.size(); i++)
    {
        file_profiling_rna << profiling_rna_mm[i][guideI][0] << "\t";
    }
    file_profiling_rna << "\n";


    //Profiling complete
    file_profiling_complete << guide << "\t";
    

    for (int i = 0; i < (len_guide); i++){ //+bulDNA
        file_profiling_complete << bp_sum_complete[i] << "\t"; //write BP columns
    }

    file_profiling_complete << "\t";
    
    int m = (len_guide + bulDNA);
    for (int i = (len_guide + bulDNA); i <  bp_sum_complete.size(); i++)
    {
        for (int j = 0; j < max_bulges; j++){
            bp_sum_complete[i] += profiling_dna_mm[m + j][guideI][0];
            //bp_sum_complete[i] += profiling_dna_mm[m + 1][guideI][0];
        }
        m += max_bulges;
    }
    m = len_guide;
    
    for (int i = (len_guide + bulDNA); i < bp_sum_complete.size(); i++)
    {
        for (int j = 0; j < max_bulges; j++){
            bp_sum_complete[i] += profiling_rna_mm[m + j][guideI][0];
            //bp_sum_complete[i] += profiling_rna_mm[m + 1][guideI][0];
        }
        m += max_bulges;
    }

    int sum_bp_complete = 0;
    for (int i = 0 ; i < len_guide; i++){ //+bulDNA
        sum_bp_complete += bp_sum_complete[i];
    }
    file_profiling_complete << ont_complete << "\t" << offt_complete << "\t" << sum_bp_complete / (double)offt_complete << "\t\t";
    for (int i = (len_guide + bulDNA); i <  bp_sum_complete.size(); i++)
        file_profiling_complete << bp_sum_complete[i] << "\t";
    
    file_profiling_complete << "\n";


    //save ext profiling
    int pos_bp_sum = len_guide + bulDNA;
    file_ext_profiling << ">" << guide << "\n";
    file_ext_profiling << "\t";
    for (int j = 0; j < len_guide; j++)
    {
        file_ext_profiling << "BP\t";
    }
    file_ext_profiling << "TARGETS\n";

    for (int mm = 0; mm <= mism; mm++)
    {
        file_ext_profiling << mm << " MISMATCHES\t";
        for (int j = 0; j < len_guide; j++)
        {
            int sum = 0;
            for (int kk = 0; kk < 4; kk++)
            {
                sum = sum + ext_profiling[guideI][mm][kk][j][0];
            }
            sum += ext_profiling_dna[guideI][mm][j][0] + ext_profiling_rna[guideI][mm][j][0];
            file_ext_profiling << sum << "\t"; //std::accumulate(extended_matrix[j][0][mm], extended_matrix[j][4][mm], 0  ) << "\t";
        }
        file_ext_profiling << bp_sum_complete[pos_bp_sum] << "\n";//profiling[len_guide + mm][guideI][0] << "\n"; // total targets
        pos_bp_sum++;
        file_ext_profiling << "NUCLEOTIDE A\t";
        for (int j = 0; j < len_guide; j++)
        {
            file_ext_profiling << ext_profiling[guideI][mm][0][j][0] << "\t";
        }
        file_ext_profiling << "\n";
        file_ext_profiling << "NUCLEOTIDE C\t";
        for (int j = 0; j < len_guide; j++)
        {
            file_ext_profiling << ext_profiling[guideI][mm][1][j][0] << "\t";
        }
        file_ext_profiling << "\n";
        file_ext_profiling << "NUCLEOTIDE G\t";
        for (int j = 0; j < len_guide; j++)
        {
            file_ext_profiling << ext_profiling[guideI][mm][2][j][0] << "\t";
        }
        file_ext_profiling << "\n";
        file_ext_profiling << "NUCLEOTIDE T\t";
        for (int j = 0; j < len_guide; j++)
        {
            file_ext_profiling << ext_profiling[guideI][mm][3][j][0] << "\t";
        }
        //file_ext_profiling << "\n\n";
        file_ext_profiling << "\n";
        file_ext_profiling << "Bulge DNA\t";
        for (int j = 0; j < len_guide; j++){
            file_ext_profiling << ext_profiling_dna[guideI][mm][j][0] << "\t";
        }

        file_ext_profiling << "\n";
        file_ext_profiling << "Bulge RNA\t";
        for (int j = 0; j < len_guide; j++){
            file_ext_profiling << ext_profiling_rna[guideI][mm][j][0] << "\t";
        }
        file_ext_profiling << "\n\n";
        
    }

}

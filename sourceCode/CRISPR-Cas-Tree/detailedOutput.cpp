#include <bitset>
#include "convert.h"
#include "detailedOutput.h"
#include "timer.h"
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
                        std::vector<std::string> &vec_guides, std::vector<std::string> &vec_targets){
    int thr = omp_get_thread_num();
    int mm_inside_string = 0;
    std::vector<std::pair<std::bitset<4>,int>> save_mm_pos; //array che salva i pair (base mm in bit, pos)
    //std::string guide_correct_dir;
    //std::string target_correct_dir;

    
    
    for (int i = 0; i < target_bit.size(); i++){ 
        //guide_correct_dir += convertBitsetToChar(guide_bit[i]);
        //target_correct_dir += convertBitsetToChar(target_bit[i]);
        

        if (i < len_guide) //check for mismatches only in guide char, not pam 
            //std::cout << "guide char: " << convert(guide_bit[i]) << ", target char: " << convert(target_bit[i]) << std::endl;
            if ((guide_bit[i] & target_bit[i]) == 0){ //mismatch case
                //target_correct_dir[i] = target_correct_dir[i] + 32 ;
                mm_inside_string++;
                //#pragma omp atomic update
                profiling[i][guideI][thr]++;
                
                save_mm_pos.push_back(std::make_pair(target_bit[i], i));
                
            }
        
    }
    if (mm_inside_string > mm) 
        std::cerr << "Il numero di mm è maggiore di quello inserito";
    //std::cout << "END TARGET" << std::endl;
    
    //std::cout << "\nguideX " << guide_correct_dir << std:: endl;
    //std::cout << "target " << target_correct_dir << std:: endl;
    //std::cout << "mm inside string: " << mm_inside_string << std::endl;
    
    //std::cout << "\nguideX " << guide_correct_dir << std:: endl;
    //std::cout << "target " << target_correct_dir << std:: endl;
    //std::cout << "mm inside string: " << mm_inside_string <<"; guideI: " << guideI << std::endl;
    //std::cout << "mm inside string: " << mm_inside_string << ", guideI: " << guideI << std::endl;
    //#pragma omp atomic update
    profiling[len_guide + mm_inside_string][guideI][thr]++; 
        //extended matrix with save_mm_pos
    for (int i = 0; i < save_mm_pos.size(); i++){
        //std::cout << "Il char mm è: " << convert(save_mm_pos[i].first) << std::endl;
        //std::cout << save_mm_pos[i].second << " " << save_mm_pos[i].second << " " << save_mm_pos[i].second << " "  << save_mm_pos[i].second << std::endl;
        //std::cout << "GuideI: "  << guideI << ", mm inside strimg: " << mm_inside_string << std::endl;
        if (save_mm_pos[i].first[0]){
            //#pragma omp atomic update
            ext_profiling[guideI][mm_inside_string][0][save_mm_pos[i].second][thr]++; //a
        }
        if (save_mm_pos[i].first[1]){
            //#pragma omp atomic update
            ext_profiling[guideI][mm_inside_string][1][save_mm_pos[i].second][thr]++; //c
        }
        if (save_mm_pos[i].first[2]){
            //#pragma omp atomic update
            ext_profiling[guideI][mm_inside_string][2][save_mm_pos[i].second][thr]++; //g
        }
        if (save_mm_pos[i].first[3]){
            //#pragma omp atomic update
            ext_profiling[guideI][mm_inside_string][3][save_mm_pos[i].second][thr]++; //t
        }
    }
    //std::cout << "guideX " << guide_correct_dir << std:: endl;
    //std::cout << "target " << target_correct_dir << std:: endl;
    
    vec_targets.push_back(t);
    vec_guides.push_back(g);
    
    
}


void detailedOutputFastBulgeDNA(int guideI, std::vector<std::bitset<4>> guide_bit, std::vector<std::bitset<4>> target_bit, std::string &g, std::string &t, 
                        int bulType, int mm, int len_guide, int bD, int bulDNA, std::vector<std::vector<std::vector<int>>> &profiling_dna, 
                        std::vector<std::vector<std::vector<int>>> &profiling_dna_mm,
                        std::vector<std::string> &vec_guides, std::vector<std::string> &vec_targets){

    int thr = omp_get_thread_num();
    //std::string guide_correct_dir;
    //std::string target_correct_dir;
    int mm_inside_string = 0;
    int bul_inside_string = 0;
    for (int i = 0; i < guide_bit.size(); i++){
        //guide_correct_dir += convertBitsetToChar(guide_bit[i]);
        //target_correct_dir += convertBitsetToChar(target_bit[i]);

        if (i < len_guide + bulDNA - bD)
            if ((guide_bit[i] & target_bit[i]) == 0){ //mismatch case
                if (guide_bit[i].none()){
                    //#pragma omp atomic update
                    profiling_dna[i][guideI][thr]++;
                    bul_inside_string++;
                    continue;
                }
                //target_correct_dir[i] = target_correct_dir[i] + 32 ;
                mm_inside_string++;
                //#pragma omp atomic update
                profiling_dna_mm[i][guideI][thr]++;
                
                
            }

    }
    //std::cout << "\nguideX " << guide_correct_dir <<  std::endl;
    //std::cout << "target " << target_correct_dir <<  std::endl;
    //std::cout << "mm inside string: " << mm_inside_string << std::endl;
    //#pragma omp atomic update
    profiling_dna_mm[len_guide + bulDNA + mm_inside_string * 2 + (bul_inside_string - 1)][guideI][thr]++; 

    
    vec_guides.push_back(g);
    vec_targets.push_back(t);
}


void detailedOutputFastBulgeRNA(int guideI, std::vector<std::bitset<4>> guide_bit, std::vector<std::bitset<4>> target_bit, std::string &g, std::string &t,
                        int bulType, int mm, int len_guide, int bD, int bulDNA, std::vector<std::vector<std::vector<int>>> &profiling_rna, 
                        std::vector<std::vector<std::vector<int>>> &profiling_rna_mm,
                        std::vector<std::string> &vec_guides, std::vector<std::string> &vec_targets){
    
    int thr = omp_get_thread_num();
    //std::cout << "thr" << thr << std::endl;
    //std::string guide_correct_dir;
    //std::string target_correct_dir;
    int mm_inside_string = 0;
    int bul_inside_string = 0;
    for (int i = 0; i < guide_bit.size(); i++){
        //guide_correct_dir += convertBitsetToChar(guide_bit[i]);
        //target_correct_dir += convertBitsetToChar(target_bit[i]);

        if (i < len_guide)
            if ((guide_bit[i] & target_bit[i]) == 0 ){ //mismatch case
                if (target_bit[i].none()){
                    //#pragma omp atomic update
                    profiling_rna[i][guideI][thr]++;
                    bul_inside_string++;
                    continue;
                }
                //target_correct_dir[i] = target_correct_dir[i] + 32 ;
                mm_inside_string++;
                //#pragma omp atomic update
                profiling_rna_mm[i][guideI][thr]++;
                
                
            }

    }
    //#pragma omp atomic update
    profiling_rna_mm[len_guide  + mm_inside_string * 2 + (bul_inside_string - 1)][guideI][thr]++; 

    //cout << "\nguideX " << guide_correct_dir <<  endl;
    //cout << "target " << target_correct_dir <<  endl;
    vec_guides.push_back(g);
    vec_targets.push_back(t);

}
void saveProfileGuide(std::string guide, int guideI, int mism, int len_guide, int bulDNA, std::vector<std::vector<std::vector<int>>> &profiling, 
                      std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling, 
                      std::vector<std::vector<std::vector<int>>> &profiling_dna, std::vector<std::vector<std::vector<int>>> &profiling_dna_mm,
                      std::vector<std::vector<std::vector<int>>> &profiling_rna, std::vector<std::vector<std::vector<int>>> &profiling_rna_mm,
                      std::ofstream &fileprofiling, std::ofstream &file_ext_profiling,
                      std::ofstream &file_profiling_dna, std::ofstream &file_profiling_rna, int num_thr){
    //sum the layer of the matrix
    for (int i = 0; i < profiling.size(); i++){
        for (int j = 1; j < num_thr; j ++){
            profiling[i][guideI][0] += profiling[i][guideI][j];
            profiling_rna[i][guideI][0] += profiling_rna[i][guideI][j];
            profiling_rna_mm[i][guideI][0] += profiling_rna_mm[i][guideI][j];
        }
    }

    for (int i = 0; i < profiling_dna_mm.size(); i++){
        for (int j = 1; j < num_thr; j ++){
            profiling_dna[i][guideI][0] += profiling_dna[i][guideI][j];
            profiling_dna_mm[i][guideI][0] += profiling_dna_mm[i][guideI][j];
        }
    }

    for (int mm = 0; mm < mism; mm++){
        for (int nuc = 0; nuc < 4; nuc++)
            for(int base = 0; base < len_guide; base++)
                for (int j = 1; j < num_thr; j++)
                    ext_profiling[guideI][mm][nuc][base][0] += ext_profiling[guideI][mm][nuc][base][j]; 
    }

    //save profiling
    fileprofiling << guide << "\t";
    for (int i = 0; i < len_guide; i++){
        fileprofiling << profiling[i][guideI][0] << "\t"; //write BP columns
    }
    fileprofiling << "\t";
    fileprofiling << profiling[len_guide][guideI][0] << "\t";//write column of ONT, take data from col with number of 0mms 
    int sum = 0;
    for (int i = (len_guide + 1); i < profiling.size(); i++){ //was 23      //sum the column 1MM,2MM...
        sum += profiling[i][guideI][0];
    }
    
    fileprofiling << sum << "\t" << profiling[len_guide][guideI][0]/(double)sum << "\t\t"; //write OFFT and ONT/OFFT
    for (int i = len_guide; i < profiling.size(); i++){     //write 0MM, 1MM etc
        fileprofiling << profiling[i][guideI][0] << "\t";
    }
    fileprofiling << "\n";

    //save ext profiling
    file_ext_profiling << ">" << guide << "\n";
    file_ext_profiling << "\t\t";
        for (int j = 0; j < len_guide; j++){
            file_ext_profiling << "BP\t";
        }
        file_ext_profiling << "TARGETS\n";
    
    for (int mm = 0; mm <= mism; mm++ ){
        file_ext_profiling << mm << "\tMISMATCHES\t";
        for (int j = 0; j < len_guide; j++){
            int sum = 0;
            for (int kk = 0; kk < 4 ; kk++){
                sum = sum + ext_profiling[guideI][mm][kk][j][0];
            }
            file_ext_profiling << sum << "\t";//std::accumulate(extended_matrix[j][0][mm], extended_matrix[j][4][mm], 0  ) << "\t";
        }
        file_ext_profiling << profiling[len_guide + mm][guideI][0] << "\n";
        file_ext_profiling << "NUCLEOTIDE\tA\t";
        for (int j = 0; j < len_guide; j++){
            file_ext_profiling << ext_profiling[guideI][mm][0][j][0] << "\t";
        }
        file_ext_profiling << "\n";
        file_ext_profiling << "NUCLEOTIDE\tC\t";
        for (int j = 0; j < len_guide; j++){
            file_ext_profiling << ext_profiling[guideI][mm][1][j][0] << "\t";
        }
        file_ext_profiling << "\n";
        file_ext_profiling << "NUCLEOTIDE\tG\t";
        for (int j = 0; j < len_guide; j++){
            file_ext_profiling << ext_profiling[guideI][mm][2][j][0] << "\t";
        }
        file_ext_profiling << "\n";
        file_ext_profiling << "NUCLEOTIDE\tT\t";
        for (int j = 0; j < len_guide; j++){
            file_ext_profiling << ext_profiling[guideI][mm][3][j][0] << "\t";
        }
        file_ext_profiling << "\n\n";
    }


    //profiling dna  
    file_profiling_dna << guide << "\t";
    for (int i = 0; i < (len_guide + bulDNA); i++){
        file_profiling_dna << profiling_dna_mm[i][guideI][0] << "(" << profiling_dna[i][guideI][0] << ")" << "\t"; //write BP columns
    }
    file_profiling_dna << "\t";
    int total_ont = profiling_dna_mm[len_guide + bulDNA][guideI][0] + profiling_dna_mm[len_guide + bulDNA + 1][guideI][0]; // sum cols 0mm1 0mm2
    file_profiling_dna << total_ont << "\t"; 
    sum = 0;
    for (int i = (len_guide + bulDNA + 2); i < profiling_dna_mm.size(); i++){ //+2 perchè devo saltare le colonne 0mm1 0mm2
        sum += profiling_dna_mm[i][guideI][0];
    }
    
    file_profiling_dna << sum << "\t" << total_ont/(double)sum << "\t\t";
    for (int i = (len_guide + bulDNA); i < profiling_dna_mm.size(); i++){
        file_profiling_dna << profiling_dna_mm[i][guideI][0] << "\t";
    }
    file_profiling_dna << "\n";


    //profiling rna
    file_profiling_rna << guide << "\t";
    for (int i = 0; i < len_guide; i++){
        file_profiling_rna << profiling_rna_mm[i][guideI][0] << "(" << profiling_rna[i][guideI][0] << ")" << "\t"; //write BP columns
    }
    file_profiling_rna << "\t";
    total_ont = profiling_rna_mm[len_guide][guideI][0] + profiling_rna_mm[len_guide + 1][guideI][0];
    file_profiling_rna << total_ont << "\t"; 
    sum = 0;
    for (int i = (len_guide + 2); i < profiling_rna_mm.size(); i++){ 
        sum += profiling_rna_mm[i][guideI][0];
    }
    
    file_profiling_rna << sum << "\t" << total_ont/(double)sum << "\t\t";
    for (int i = (len_guide); i < profiling_rna_mm.size(); i++){
        file_profiling_rna << profiling_rna_mm[i][guideI][0] << "\t";
    }
    file_profiling_rna << "\n";


}
#include <bitset>
#include "convert.h"
#include "detailedOutput.h"
#include "timer.h"

#include <map>
#include <numeric> //for accumulate

/**This function provides a profiling of the results from the searchOnTST. The results are computed
 * as soon as the target is found on the TST
 * @param guideI: index of the guide 
 * @param guide: the string containing the guide
 * @param target: the string containing the target found
 * @param bulType: int value of the bulge type (0 -> " X ", <0 -> "RNA", >0 -> "DNA")
 * @param mm: number of mismatches allowed
 * @param profiling: an 23+mm x numGuides matrix for profiling.Cells 0-19 are for counting positional
 * mms, cell 20 for ONT, cell 21 for OFFT,  cell 22+ for total string count mms
 * devo farmi passare usando & anche le matrici di salvataggio finale
 * @param ext_profiling: matrix for extended profiling. The first dimension selects the guideI, the second dimension selects
 * the matrix to fill w.r.t. the number of mms in the target, the third dimension select the nucleotide (A,C,G,T) and the 
 * fourth dimension selects the position in which there's the mm
 */
void detailedOutputFast(int guideI, std::vector<std::bitset<4>> guide_bit, std::vector<std::bitset<4>> target_bit, 
                        int bulType, int mm, std::vector<std::vector<int>> &profiling, std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling,
                        std::vector<std::string> &vec_guides, std::vector<std::string> &vec_targets){
    if (bulType != 0)       //TODO add the conversion guide_bit -> guide; target_bit->target and add to result for RNA-DNA bulges
        return;
    
    int mm_inside_string = 0;
    std::vector<std::pair<std::bitset<4>,int>> save_mm_pos; //array che salva i pair (base mm in bit, pos)
    std::string guide_correct_dir;
    std::string target_correct_dir; 
    
    for (int i = 0; i < target_bit.size(); i++){ 
        guide_correct_dir += convert(guide_bit[i]);
        target_correct_dir += convert(target_bit[i]);
        if (i < 20)
        //std::cout << "guide char: " << convert(guide_bit[i]) << ", target char: " << convert(target_bit[i]) << std::endl;
        if ((guide_bit[i] & target_bit[i]) == 0){ //mismatch case
            target_correct_dir[i] = target_correct_dir[i] + 32 ;
            mm_inside_string++;
            profiling[i][guideI]++;
            
            
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
    profiling[22 + mm_inside_string][guideI]++; 
        //extended matrix with save_mm_pos
    for (int i = 0; i < save_mm_pos.size(); i++){
        //std::cout << "Il char mm è: " << convert(save_mm_pos[i].first) << std::endl;
        //std::cout << save_mm_pos[i].second << " " << save_mm_pos[i].second << " " << save_mm_pos[i].second << " "  << save_mm_pos[i].second << std::endl;
        //std::cout << "GuideI: "  << guideI << ", mm inside strimg: " << mm_inside_string << std::endl;
        if (save_mm_pos[i].first[0])
            ext_profiling[guideI][mm_inside_string][0][save_mm_pos[i].second]++; //a
        if (save_mm_pos[i].first[1])
            ext_profiling[guideI][mm_inside_string][1][save_mm_pos[i].second]++; //c
        if (save_mm_pos[i].first[2])
            ext_profiling[guideI][mm_inside_string][2][save_mm_pos[i].second]++; //g
        if (save_mm_pos[i].first[3])
            ext_profiling[guideI][mm_inside_string][3][save_mm_pos[i].second]++; //t
    }
    //std::cout << "guideX " << guide_correct_dir << std:: endl;
    //std::cout << "target " << target_correct_dir << std:: endl;
    
    vec_targets.push_back(target_correct_dir);
    vec_guides.push_back(guide_correct_dir);
    
    
}

void saveProfileGuide(std::string guide, int guideI, int mism, std::vector<std::vector<int>> &profiling, 
                      std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling, 
                      std::ofstream &fileprofiling, std::ofstream &file_ext_profiling){
    
    //save profiling
    fileprofiling << guide << "\t";
    for (int i = 0; i < 20; i++){
        fileprofiling << profiling[i][guideI] << "\t";
    }
    fileprofiling << "\t";
    fileprofiling << profiling[22][guideI] << "\t";
    int sum = 0;
    for (int i = 23; i < profiling.size(); i++){
        sum += profiling[i][guideI];
    }
    
    fileprofiling << sum << "\t" << profiling[22][guideI]/(double)sum << "\t\t";
    for (int i = 22; i < profiling.size(); i++){
        fileprofiling << profiling[i][guideI] << "\t";
    }
    fileprofiling << "\n";

    //save ext profiling
    file_ext_profiling << ">" << guide << "\n";
    file_ext_profiling << "\t\t";
        for (int j = 0; j < 20; j++){
            file_ext_profiling << "BP\t";
        }
        file_ext_profiling << "TARGETS\n";
    
    for (int mm = 0; mm <= mism; mm++ ){
        file_ext_profiling << mm << "\tMISMATCHES\t";
        for (int j = 0; j < 20; j++){
            int sum = 0;
            for (int kk = 0; kk < 4 ; kk++){
                sum = sum + ext_profiling[guideI][mm][kk][j];
            }
            file_ext_profiling << sum << "\t";//std::accumulate(extended_matrix[j][0][mm], extended_matrix[j][4][mm], 0  ) << "\t";
        }
        file_ext_profiling << profiling[22 + mm][guideI] << "\n";
        file_ext_profiling << "NUCLEOTIDE\tA\t";
        for (int j = 0; j < 20; j++){
            file_ext_profiling << ext_profiling[guideI][mm][0][j] << "\t";
        }
        file_ext_profiling << "\n";
        file_ext_profiling << "NUCLEOTIDE\tC\t";
        for (int j = 0; j < 20; j++){
            file_ext_profiling << ext_profiling[guideI][mm][1][j] << "\t";
        }
        file_ext_profiling << "\n";
        file_ext_profiling << "NUCLEOTIDE\tG\t";
        for (int j = 0; j < 20; j++){
            file_ext_profiling << ext_profiling[guideI][mm][2][j] << "\t";
        }
        file_ext_profiling << "\n";
        file_ext_profiling << "NUCLEOTIDE\tT\t";
        for (int j = 0; j < 20; j++){
            file_ext_profiling << ext_profiling[guideI][mm][3][j] << "\t";
        }
        file_ext_profiling << "\n\n";
    }
}
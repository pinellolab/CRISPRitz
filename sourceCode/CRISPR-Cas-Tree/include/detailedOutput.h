#pragma once

#include <iostream>
#include <fstream>
#include <vector>



void detailedOutputFast(int guideI, std::vector<std::bitset<4>> guide_bit, std::vector<std::bitset<4>> target_bit, std::string &g, std::string &t,
                        int bulType, int mm, int len_guide, std::vector<std::vector<std::vector<int>>> &profiling, std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling,
                        std::vector<std::string> &vec_guides, std::vector<std::string> &vec_targets, bool pam_at_start); 
        
void detailedOutputFastBulgeDNA(int guideI, std::vector<std::bitset<4>> guide_bit, std::vector<std::bitset<4>> target_bit, std::string &g, std::string &t,
                        int bulType, int mm, int max_bulges, int len_guide, int bD, int bulDNA, std::vector<std::vector<std::vector<int>>> &profiling_dna, 
                        std::vector<std::vector<std::vector<int>>> &profiling_dna_mm,
                        std::vector<std::string> &vec_guides, std::vector<std::string> &vec_targets, bool pam_at_start,
                        std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_dna,
                        std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling); 

void detailedOutputFastBulgeRNA(int guideI, std::vector<std::bitset<4>> guide_bit, std::vector<std::bitset<4>> target_bit, std::string &g, std::string &t,
                        int bulType, int mm, int max_bulges, int len_guide, int bD, int bulDNA, std::vector<std::vector<std::vector<int>>> &profiling_rna, 
                        std::vector<std::vector<std::vector<int>>> &profiling_rna_mm,
                        std::vector<std::string> &vec_guides, std::vector<std::string> &vec_targets, bool pam_at_start,
                        std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profilingrna,
                        std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling); 

void saveProfileGuide(std::string guide, int guideI, int mism, int max_bulges, int len_guide, int bulDNA, std::vector<std::vector<std::vector<int>>> &profiling, 
                      std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> &ext_profiling, 
                      std::vector<std::vector<std::vector<int>>> &profiling_dna, std::vector<std::vector<std::vector<int>>> &profiling_dna_mm,
                      std::vector<std::vector<std::vector<int>>> &profiling_rna, std::vector<std::vector<std::vector<int>>> &profiling_rna_mm,
                      std::ofstream &fileprofiling, std::ofstream &file_ext_profiling,
                      std::ofstream &file_profiling_dna, std::ofstream &file_profiling_rna, std::ofstream &file_profiling_complete, int num_thr, int pam_size, bool pam_at_start,
                      std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_dna,
                      std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling_rna);

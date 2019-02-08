#pragma once

#include <iostream>
#include <fstream>
#include <vector>

void detailedOutputFast(int guideI, std::vector<std::bitset<4>> guide_bit, std::vector<std::bitset<4>> target_bit, 
                        int bulType, int mm, std::vector<std::vector<int>> &profiling, std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling, 
                        std::vector<std::string> &vec_guides, std::vector<std::string> &vec_targets); 

void saveProfileGuide(std::string guide, int guideI, int mism, std::vector<std::vector<int>> &profiling,
                      std::vector<std::vector<std::vector<std::vector<int>>>> &ext_profiling, 
                      std::ofstream &fileprofiling, std::ofstream &file_ext_profiling);
#pragma once

#include <iostream>
#include <string>
#include <vector>

struct data_str{
    std::string bT, vIG, vTOG, chr;
    int ind, mm, bul;
    char dir;
};

bool compareIndices(const data_str a, const data_str b);

void resultFilter(std::vector<std::string> &bulgeType, std::vector<std::string> &vecInGuide, std::vector<std::string> &vecTargetOfGuide,
                  std::string &chrName, std::vector<int> &indices, std::vector<int> &mismatches, std::vector<int> &bulgeSize, 
                  std::vector<char> &directions);
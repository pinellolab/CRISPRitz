#include <algorithm>
#include <chrono>	//for time measurements


#include "filtering.h"

// struct data_str{
//     std::string bT, vIG, vTOG, chr;
//     int ind, mm, bul;
//     char dir;
// }

bool compareIndices(const data_str a, const data_str b){
    return a.ind < b.ind;    
}

void resultFilter(std::vector<std::string> &bulgeType, std::vector<std::string> &vecInGuide, std::vector<std::string> &vecTargetOfGuide,
                  std::string &chrName, std::vector<int> &indices, std::vector<int> &mismatches, std::vector<int> &bulgeSize, 
                  std::vector<char> &directions){
    
    std::vector<data_str> pre_filter;
    std::vector<data_str> post_filter;

    for (int i = 0; i < indices.size(); i++){
        data_str d = {bulgeType[i], vecInGuide[i] , vecTargetOfGuide[i] , chrName , indices[i]
                            , mismatches[i] ,bulgeSize[i], directions[i]  };
        if (bulgeType[i].find("X") != std::string::npos) //è X                    
            post_filter.push_back(d);
        else{
            pre_filter.push_back(d);
        }

    }

    std::sort(pre_filter.begin(), pre_filter.end(), compareIndices);

    data_str min = pre_filter[0];

    //Filter only the DNA/RNA bulges targets
    for (int i = 1; i < pre_filter.size(); i++){
        if (min.ind == pre_filter[i].ind){
            do{
                if (min.mm < pre_filter[i].mm && min.bul < pre_filter[i].bul){
                    min = pre_filter[i];
                }
                i++;
            }while(i < pre_filter.size() && min.ind == pre_filter[i].ind);
        }
        post_filter.push_back(min);
        if (i < pre_filter.size())
            min = pre_filter[i];
    }

}
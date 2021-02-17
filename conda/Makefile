#compile crispritz cpp softwares
crispritz:
	$(CXX) -std=c++11 -O3 -fopenmp -I $(BUILD_PREFIX)/include sourceCode/CRISPR-Cas-Tree/main.cpp -o buildTST
	$(CXX) -std=c++11 -O3 -fopenmp -I $(BUILD_PREFIX)/include sourceCode/CRISPR-Cas-Tree/searchOnTST.cpp sourceCode/CRISPR-Cas-Tree/detailedOutput.cpp sourceCode/CRISPR-Cas-Tree/convert.cpp -I sourceCode/CRISPR-Cas-Tree/include -o searchTST
	$(CXX) -std=c++11 -O3 -fopenmp -I $(BUILD_PREFIX)/include sourceCode/CRISPRofiler/main.cpp sourceCode/CRISPRofiler/profiling.cpp sourceCode/CRISPRofiler/guide_searching.cpp sourceCode/CRISPRofiler/pam_searching.cpp sourceCode/CRISPRofiler/pre_computation.cpp sourceCode/CRISPRofiler/reading.cpp sourceCode/CRISPRofiler/analysis.cpp -o searchBruteForce

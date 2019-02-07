#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include <ctype.h>
#include <vector>
#include <algorithm>
#include <dirent.h>
#include <unistd.h>
#include <cstring>
#include <fcntl.h>
#include <cmath>
#include <iterator>
#include <omp.h>
#include <bits/stdc++.h>
#include <bitset>

using namespace std;

void analysis();
void reading_pam();
void reading_guide();
void reading_chromosomes(char **argv);
void reverse_guide();
void generateprimenumbers();
int buildMatchingMachine(string arr[], int k);
void searchPam();
string switchSymbol(char sym);
string reverse(string pamInput);
vector<string> getProducts(string s[], int s_size);
vector<string> generatePam(string pamInput);
string reversetarget(string targetfound);
string missmatching(string targetfound,string guidapassata,int guidafound);
void checkguide(int guidaprima);
void guide_searching();
void profiler();
void read_bed();
void generate_tree(bool flag);
void interval_profiler();
void genomebitconversion();
void guidesbitconversion();
void profilersetting();
void resetter(int inizio);
void read_special(string path);
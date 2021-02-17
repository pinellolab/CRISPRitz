#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <queue>
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
//#include <bits/stdc++.h>
#include <bitset>
#include <algorithm>
#include <iomanip>
#include <boost/dynamic_bitset.hpp>

using namespace std;

//dimensions for automata
#define MAXS 100000
#define MAXC 93
#define MAXW 100000

//MIN and MAX definition
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

//function definition
void analysis();
void guide_searching();
string missmatching(string targetfound, int index, int guidafound, int reverse);
string reversetarget(string targetfound);
int main(int argc, char **argv);
void buildMachine();
vector<string> generatePam(string pamInput);
vector<string> getProducts(string s[], int s_size);
string reverse(string pamInput);
void searchPam();
string switchSymbol_novar(char sym);
string switchSymbol_var(char sym);
void genomebitconversion();
void guidesbitconversion();
void profilersetting();
void profiler();
void reading_chromosomes(char **argv);
void reading_guide();
void reading_pam();
void pamGeneration();

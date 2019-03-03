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
#include "crispritz.h"

using namespace std;

extern int genlen;
extern double totaltimepam,totaltimeguide;
extern vector<int> pamindices,pamindicesreverse;
extern vector<bitset<4>> genomebit;
extern string genome;

void analysis()
{
   //variable reset
   pamindices.clear();
   pamindicesreverse.clear();
   //pamindices.resize(genlen);
   //pamindicesreverse.resize(genlen);

   //start pam searching
   searchPam();

   //start guides searching, executed number of guides times
   guide_searching();

   //clear the genome string for the next genome analysis
   genome.clear();
   genomebit.clear();
}
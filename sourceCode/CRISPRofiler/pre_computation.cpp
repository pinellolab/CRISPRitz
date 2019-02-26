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

extern int i, j, genlen, guidelen, pamlimit, totalguides, threads, guidelencorrected, inputmissmatch;
extern vector<string> guides;
extern string guide, genome;
extern vector<vector<int>> guideprofiling;                  //vettore contenente il profilo
extern vector<vector<vector<vector<int>>>> matrixprofiling; //vettore contenente il profilo extended
extern vector<bitset<4>> genomebit;                         //genoma in bit
extern vector<vector<bitset<4>>> guidesbit;                 //guide in bit
extern vector<vector<bitset<4>>> reverseguidesbit;          //reverse guide in bit

void genomebitconversion() //converto il genoma dal fasta alla versione bit
{
   genomebit.clear();
   genomebit.resize(genlen);

#pragma omp parallel for num_threads(threads) schedule(static)
   for (i = 0; i < genlen; i++)
   {
      if (genome[i] == 'A')
      {
         genomebit[i] = bitset<4>(string("0001"));
      }
      else if (genome[i] == 'C')
      {
         genomebit[i] = bitset<4>(string("0010"));
      }
      else if (genome[i] == 'G')
      {
         genomebit[i] = bitset<4>(string("0100"));
      }
      else if (genome[i] == 'T')
      {
         genomebit[i] = bitset<4>(string("1000"));
      }
      else if (genome[i] == 'N')
      {
         genomebit[i] = bitset<4>(string("0000"));
      }
      else if (genome[i] == 'R')
      {
         genomebit[i] = bitset<4>(string("0101"));
      }
      else if (genome[i] == 'Y')
      {
         genomebit[i] = bitset<4>(string("1010"));
      }
      else if (genome[i] == 'S')
      {
         genomebit[i] = bitset<4>(string("0110"));
      }
      else if (genome[i] == 'W')
      {
         genomebit[i] = bitset<4>(string("1001"));
      }
      else if (genome[i] == 'K')
      {
         genomebit[i] = bitset<4>(string("1100"));
      }
      else if (genome[i] == 'M')
      {
         genomebit[i] = bitset<4>(string("0011"));
      }
      else if (genome[i] == 'B')
      {
         genomebit[i] = bitset<4>(string("1110"));
      }
      else if (genome[i] == 'D')
      {
         genomebit[i] = bitset<4>(string("1101"));
      }
      else if (genome[i] == 'H')
      {
         genomebit[i] = bitset<4>(string("1011"));
      }
      else if (genome[i] == 'V')
      {
         genomebit[i] = bitset<4>(string("0111"));
      }
   }
}

void guidesbitconversion() //converto le guide dal file.txt alla versione bit
{
   guidesbit.resize(totalguides);
   reverseguidesbit.resize(totalguides);
   guidelen = guides[0].size();
   guidelencorrected = guidelen - pamlimit;

#pragma omp parallel for num_threads(threads) schedule(static) private(j)
   for (i = 0; i < totalguides; i++)
   {
      guidesbit[i].resize(guidelen);
      reverseguidesbit[i].resize(guidelen);

      for (j = 0; j < guidelencorrected; j++)
      {
         if (guides[i][j] == 'A')
         {
            guidesbit[i][j] = bitset<4>(string("0001"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("1000"));
         }
         else if (guides[i][j] == 'C')
         {
            guidesbit[i][j] = bitset<4>(string("0010"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("0100"));
         }
         else if (guides[i][j] == 'G')
         {

            guidesbit[i][j] = bitset<4>(string("0100"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("0010"));
         }
         else if (guides[i][j] == 'T')
         {
            guidesbit[i][j] = bitset<4>(string("1000"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("0001"));
         }
         else if (guides[i][j] == 'N')
         {
            guidesbit[i][j] = bitset<4>(string("1111"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("1111"));
         }
         else if (guides[i][j] == 'R')
         {
            guidesbit[i][j] = bitset<4>(string("0101"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("1010"));
         }
         else if (guides[i][j] == 'Y')
         {
            guidesbit[i][j] = bitset<4>(string("1010"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("0101"));
         }
         else if (guides[i][j] == 'S')
         {
            guidesbit[i][j] = bitset<4>(string("0110"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("1001"));
         }
         else if (guides[i][j] == 'W')
         {
            guidesbit[i][j] = bitset<4>(string("1001"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("0110"));
         }
         else if (guides[i][j] == 'K')
         {
            guidesbit[i][j] = bitset<4>(string("1100"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("0011"));
         }
         else if (guides[i][j] == 'M')
         {
            guidesbit[i][j] = bitset<4>(string("0011"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("1100"));
         }
         else if (guides[i][j] == 'B')
         {
            guidesbit[i][j] = bitset<4>(string("1110"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("0111"));
         }
         else if (guides[i][j] == 'D')
         {
            guidesbit[i][j] = bitset<4>(string("1101"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("1011"));
         }
         else if (guides[i][j] == 'H')
         {
            guidesbit[i][j] = bitset<4>(string("1011"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("1101"));
         }
         else if (guides[i][j] == 'V')
         {
            guidesbit[i][j] = bitset<4>(string("0111"));
            reverseguidesbit[i][guidelen - 1 - j] = bitset<4>(string("1110"));
         }
      }
   }
}

void profilersetting() //alloco spazio per matrici dei profili
{
   guideprofiling.resize(totalguides);
   matrixprofiling.resize(totalguides);

   for (i = 0; i < totalguides; i++)
   {
      guideprofiling[i].resize(guidelen - pamlimit + 3 + inputmissmatch + 1);
      matrixprofiling[i].resize(inputmissmatch + 1);

      for (int kk = 0; kk <= inputmissmatch; kk++)
      {
         matrixprofiling[i][kk].resize(guidelen - pamlimit);

         for (int jj = 0; jj < guidelen - pamlimit; jj++)
         {
            matrixprofiling[i][kk][jj].resize(5);
         }
      }
   }
}

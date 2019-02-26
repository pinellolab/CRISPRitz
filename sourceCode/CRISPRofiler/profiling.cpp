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

extern int i, j, guidelen, pamlimit, totalguides;
extern vector<string> guides;
extern vector<vector<int>> guideprofiling;
extern vector<vector<vector<vector<int>>>> matrixprofiling;
extern string writeprofile, writeextensiveprofile;
extern int guidelen,inputmissmatch,pamlimit;

void profiler()
{
   //writeprofile = "GUIDE\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\t\tONT\tOFFT\tOFFT/MM\t\t0MM\t1MM\t2MM\t3MM\t4MM\t5MM\t6MM\t7MM\t8MM\t9MM\t10MM\n";
   writeprofile ="GUIDE\t";
   for (int ff=0;ff<(guidelen-pamlimit);ff++)
   {
      writeprofile += "BP\t";
   }

   writeprofile += "\tONT\tOFFT\tOFFT/MM\t\t";

   for (int ff=0;ff<=inputmissmatch;ff++)
   {
      writeprofile += to_string(ff)+"MM\t";
   }

   writeprofile += "\n";
   

   for (i = 0; i < totalguides; i++)
   {
      writeprofile += guides[i]; //scrivo la guida iesima nel profilo
      writeprofile += "\t";

      for (j = 0; j < (guidelen - pamlimit); j++)
      {
         writeprofile += to_string(guideprofiling[i][j]); //scrivo tutti i valori delle basi (mismatches generati da ogni base)
         writeprofile += "\t";
      }

      writeprofile += "\t";
      writeprofile += to_string(guideprofiling[i][guidelen-pamlimit]); //on-target (off-target con 0 mismatch)

      writeprofile += "\t";
      writeprofile += to_string(guideprofiling[i][(guidelen-pamlimit)+1]); //off-target
      writeprofile += "\t";

      double offtarget = guideprofiling[i][(guidelen-pamlimit)+1]; // offtarget
      double missmatch = guideprofiling[i][(guidelen-pamlimit)+2]; // offtarget/missmatch
      double meanmissmatch = missmatch / offtarget;

      writeprofile += to_string(meanmissmatch);
      writeprofile += "\t\t";

      for (j = (guidelen - pamlimit + 3); j <= (guidelen - pamlimit + 3 + inputmissmatch); j++)
      {
         writeprofile += to_string(guideprofiling[i][j]); //scrivo tutti i risultati delle colonne con threshold di mismatches
         writeprofile += "\t";
      }

      writeprofile += "\n";

      writeextensiveprofile += ">" + guides[i];
      writeextensiveprofile += "\n";
      writeextensiveprofile += "\t";
      //writeextensiveprofile += "\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tTARGETS\n";

      for (int ff=0;ff<(guidelen-pamlimit);ff++)
      {
         writeextensiveprofile += "BP\t";
      }

      writeextensiveprofile += "TARGETS";
      writeextensiveprofile += "\n";

      for (int hh = 0; hh <=inputmissmatch; hh++)
      {
         writeextensiveprofile += to_string(hh) + " MISMATCHES"; //scrivo le righe 0 MISMATCH, 1 MISMATCH, ECC

         for (int dd = 0; dd < (guidelen - pamlimit); dd++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(matrixprofiling[i][hh][dd][4]); // totale mismatches per base (in quella threshold di mismatches)
         }

         writeextensiveprofile += "\t";
         writeextensiveprofile += to_string(guideprofiling[i][(guidelen-pamlimit)+3 + hh]); //totale off-targets in quella threshold di mismatches

         writeextensiveprofile += "\n";

         writeextensiveprofile += "NUCLEOTIDE A";

         for (int ff = 0; ff < (guidelen - pamlimit); ff++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(matrixprofiling[i][hh][ff][0]); //totale mismatches per BP con nucleotide indicato
         }

         writeextensiveprofile += "\n";
         writeextensiveprofile += "NUCLEOTIDE C";

         for (int ff = 0; ff < (guidelen - pamlimit); ff++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(matrixprofiling[i][hh][ff][1]); //totale mismatches per BP con nucleotide indicato
         }

         writeextensiveprofile += "\n";
         writeextensiveprofile += "NUCLEOTIDE G";

         for (int ff = 0; ff < (guidelen - pamlimit); ff++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(matrixprofiling[i][hh][ff][2]); //totale mismatches per BP con nucleotide indicato
         }

         writeextensiveprofile += "\n";
         writeextensiveprofile += "NUCLEOTIDE T";

         for (int ff = 0; ff < (guidelen - pamlimit); ff++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(matrixprofiling[i][hh][ff][3]); //totale mismatches per BP con nucleotide indicato
         }

         writeextensiveprofile += "\n";
         writeextensiveprofile += "\n";
      }
   }
}
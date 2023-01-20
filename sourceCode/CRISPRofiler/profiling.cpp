#include "include/crispritz.h"

using namespace std;

extern int guidelen, pamlimit, totalguides, inputmissmatch;
extern vector<string> guides;
extern vector<vector<int>> guideprofiling;
extern vector<vector<vector<vector<int>>>> matrixprofiling;
extern string writeprofile, writeextensiveprofile;

void profiler()
{
   // writeprofile = "GUIDE\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\t\tONT\tOFFT\tOFFT/MM\t\t0MM\t1MM\t2MM\t3MM\t4MM\t5MM\t6MM\t7MM\t8MM\t9MM\t10MM\n";
   writeprofile = "GUIDE\t";
   for (int ff = 0; ff < guidelen; ff++)
   {
      writeprofile += "BP\t";
   }

   writeprofile += "\tONT\tOFFT\tOFFT/MM\t\t";

   for (int ff = 0; ff <= inputmissmatch; ff++)
   {
      writeprofile += to_string(ff) + "MM\t";
   }

   writeprofile += "\n";

   for (int i = 0; i < guides.size(); i++)
   {
      writeprofile += guides[i]; // scrivo la guida iesima nel profilo
      writeprofile += "\t";

      for (int j = 0; j < guidelen; j++)
      {
         writeprofile += to_string(guideprofiling[i][j]); // scrivo tutti i valori delle basi (mismatches generati da ogni base)
         writeprofile += "\t";
      }

      writeprofile += "\t";
      writeprofile += to_string(guideprofiling[i][guidelen]); // on-target (off-target con 0 mismatch)

      writeprofile += "\t";
      writeprofile += to_string(guideprofiling[i][(guidelen) + 1]); // off-target
      writeprofile += "\t";

      double offtarget = guideprofiling[i][(guidelen) + 1]; // offtarget
      double missmatch = guideprofiling[i][(guidelen) + 2]; // offtarget/missmatch
      double meanmissmatch = missmatch / offtarget;

      writeprofile += to_string(meanmissmatch);
      writeprofile += "\t\t";

      for (int j = guidelen + 3; j <= (guidelen + 3 + inputmissmatch); j++)
      {
         writeprofile += to_string(guideprofiling[i][j]); // scrivo tutti i risultati delle colonne con threshold di mismatches
         writeprofile += "\t";
      }

      writeprofile += "\n";

      writeextensiveprofile += ">" + guides[i];
      writeextensiveprofile += "\n";
      writeextensiveprofile += "\t";
      // writeextensiveprofile += "\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tTARGETS\n";

      for (int ff = 0; ff < (guidelen); ff++)
      {
         writeextensiveprofile += "BP\t";
      }

      writeextensiveprofile += "TARGETS";
      writeextensiveprofile += "\n";

      for (int current_mm = 0; current_mm <= inputmissmatch; current_mm++)
      {
         writeextensiveprofile += to_string(current_mm) + " MISMATCHES"; // scrivo le righe 0 MISMATCH, 1 MISMATCH, ECC

         for (int guide_nt = 0; guide_nt < (guidelen); guide_nt++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(matrixprofiling[i][current_mm][guide_nt][4]); // totale mismatches per base (in quella threshold di mismatches)
         }

         writeextensiveprofile += "\t";
         writeextensiveprofile += to_string(guideprofiling[i][guidelen + 3 + current_mm]); // totale off-targets in quella threshold di mismatches

         writeextensiveprofile += "\n";

         writeextensiveprofile += "NUCLEOTIDE A";

         for (int ff = 0; ff < (guidelen); ff++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(matrixprofiling[i][current_mm][ff][0]); // totale mismatches per BP con nucleotide indicato
         }

         writeextensiveprofile += "\n";
         writeextensiveprofile += "NUCLEOTIDE C";

         for (int ff = 0; ff < (guidelen); ff++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(matrixprofiling[i][current_mm][ff][1]); // totale mismatches per BP con nucleotide indicato
         }

         writeextensiveprofile += "\n";
         writeextensiveprofile += "NUCLEOTIDE G";

         for (int ff = 0; ff < (guidelen); ff++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(matrixprofiling[i][current_mm][ff][2]); // totale mismatches per BP con nucleotide indicato
         }

         writeextensiveprofile += "\n";
         writeextensiveprofile += "NUCLEOTIDE T";

         for (int ff = 0; ff < (guidelen); ff++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(matrixprofiling[i][current_mm][ff][3]); // totale mismatches per BP con nucleotide indicato
         }

         writeextensiveprofile += "\n";
         writeextensiveprofile += "Bulge DNA";

         for (int ff = 0; ff < (guidelen); ff++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(0); // totale mismatches per BP con nucleotide indicato
         }

         writeextensiveprofile += "\n";
         writeextensiveprofile += "Bulge RNA";

         for (int ff = 0; ff < (guidelen); ff++)
         {
            writeextensiveprofile += "\t";
            writeextensiveprofile += to_string(0); // totale mismatches per BP con nucleotide indicato
         }

         writeextensiveprofile += "\n";
         writeextensiveprofile += "\n";
      }
   }
}

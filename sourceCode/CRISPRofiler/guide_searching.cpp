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

extern ofstream results;
extern int i, j, guidelen, currentmissmatch, pamlimit, guidecount, totalguides, threads, countmissmatchpos, guidelencorrected;
extern double totaltimeguidesearch;
extern vector<int> missmatchthrestotal, pamindices, pamindicesreverse;
extern vector<string> guides, guidestringlist;
extern string chrnames, guide, genome, totalbuf;
extern char nowriting, noprofile;

//profile
extern vector<vector<int>> guideprofiling;
extern vector<vector<vector<vector<int>>>> matrixprofiling;

//genome
extern vector<bitset<4>> genomebit, targetfoundbit;
extern vector<vector<bitset<4>>> guidesbit;
extern vector<vector<bitset<4>>> reverseguidesbit;

string reversetarget(string targetfound) //reverse del target per poterlo confrontare con la guida
{
   targetfoundbit.clear();
   targetfoundbit.resize(guidelen);

   for (j = 0; j < guidelen; j++)
   {
      if (targetfound[j] == 'A')
      {
         targetfound[j] = 'T';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("1000"));
      }
      else if (targetfound[j] == 'T')
      {
         targetfound[j] = 'A';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("0001"));
      }
      else if (targetfound[j] == 'G')
      {
         targetfound[j] = 'C';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("0010"));
      }
      else if (targetfound[j] == 'C')
      {
         targetfound[j] = 'G';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("0100"));
      }
      else if (targetfound[j] == 'R')
      {
         targetfound[j] = 'Y';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("1010"));
      }
      else if (targetfound[j] == 'Y')
      {
         targetfound[j] = 'R';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("0101"));
      }
      else if (targetfound[j] == 'S')
      {
         targetfound[j] = 'S';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("0110"));
      }
      else if (targetfound[j] == 'W')
      {
         targetfound[j] = 'W';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("1001"));
      }
      else if (targetfound[j] == 'M')
      {
         targetfound[j] = 'K';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("1100"));
      }
      else if (targetfound[j] == 'K')
      {
         targetfound[j] = 'M';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("0011"));
      }
      else if (targetfound[j] == 'H')
      {
         targetfound[j] = 'D';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("1101"));
      }
      else if (targetfound[j] == 'D')
      {
         targetfound[j] = 'H';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("1011"));
      }
      else if (targetfound[j] == 'B')
      {
         targetfound[j] = 'V';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("0111"));
      }
      else if (targetfound[j] == 'V')
      {
         targetfound[j] = 'B';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("1110"));
      }
      else if (targetfound[j] == 'N')
      {
         targetfound[j] = 'N';
         targetfoundbit[guidelen - 1 - j] = bitset<4>(string("1111"));
      }
   }

   reverse(targetfound.begin(), targetfound.end());
   return targetfound;
}

string missmatching(string targetfound, int index, int guidafound, int reverse) //conteggio dei mismatches e profilazione se richiesta dall'utente
{
   if (noprofile == 'p') //se richiesto profilo entro
   {
      //assegnazione valori ai contatori temporanei
      vector<char> contamissmatch;
      vector<vector<float>> contabasi;
      countmissmatchpos = 0;
      contamissmatch.resize(guidelen - pamlimit);
      contabasi.resize(guidelen - pamlimit);

      if (reverse == 0)
      {
         //inizio ricerca sulla guida
         for (j = 0; j < guidelencorrected; j++)
         {
            contabasi[j].resize(4); //contatore per motiflogo

            if ((genomebit[index + j] & guidesbit[guidafound][j]) == 0)
            {
               //catena di elseif per verificare quale base genera MM
               if (targetfound[j] == 'A')
                  contabasi[j][0]++;
               else if (targetfound[j] == 'T')
                  contabasi[j][3]++;
               else if (targetfound[j] == 'G')
                  contabasi[j][2]++;
               else if (targetfound[j] == 'C')
                  contabasi[j][1]++;
               else if (targetfound[j] == 'R')
               {
                  contabasi[j][0] = contabasi[j][0] + 0.5;
                  contabasi[j][2] = contabasi[j][2] + 0.5;
               }
               else if (targetfound[j] == 'Y')
               {
                  contabasi[j][1] = contabasi[j][1] + 0.5;
                  contabasi[j][3] = contabasi[j][3] + 0.5;
               }
               else if (targetfound[j] == 'S')
               {
                  contabasi[j][2] = contabasi[j][2] + 0.5;
                  contabasi[j][1] = contabasi[j][1] + 0.5;
               }
               else if (targetfound[j] == 'W')
               {
                  contabasi[j][0] = contabasi[j][0] + 0.5;
                  contabasi[j][3] = contabasi[j][3] + 0.5;
               }
               else if (targetfound[j] == 'K')
               {
                  contabasi[j][2] = contabasi[j][2] + 0.5;
                  contabasi[j][3] = contabasi[j][3] + 0.5;
               }
               else if (targetfound[j] == 'M')
               {
                  contabasi[j][0] = contabasi[j][0] + 0.5;
                  contabasi[j][1] = contabasi[j][1] + 0.5;
               }
               else if (targetfound[j] == 'B')
               {
                  contabasi[j][1] = contabasi[j][1] + (1 / 3);
                  contabasi[j][2] = contabasi[j][2] + (1 / 3);
                  contabasi[j][3] = contabasi[j][3] + (1 / 3);
               }
               else if (targetfound[j] == 'D')
               {
                  contabasi[j][0] = contabasi[j][0] + (1 / 3);
                  contabasi[j][2] = contabasi[j][2] + (1 / 3);
                  contabasi[j][3] = contabasi[j][3] + (1 / 3);
               }
               else if (targetfound[j] == 'H')
               {
                  contabasi[j][0] = contabasi[j][0] + (1 / 3);
                  contabasi[j][1] = contabasi[j][1] + (1 / 3);
                  contabasi[j][3] = contabasi[j][3] + (1 / 3);
               }
               else if (targetfound[j] == 'V')
               {
                  contabasi[j][0] = contabasi[j][0] + (1 / 3);
                  contabasi[j][1] = contabasi[j][1] + (1 / 3);
                  contabasi[j][2] = contabasi[j][2] + (1 / 3);
               }

               targetfound[j] = tolower(targetfound[j]); //metto lettera minuscola su guida da scrivere in output
               countmissmatchpos++;                      //conto MM totali della guida
               contamissmatch[j] = 5;                    //segno posizione MM
               guideprofiling[guidafound][j]++;          //segno un MM in una base della guida nel profilo
            }
         }
      }
      else
      {
         for (j = 0; j < guidelencorrected; j++)
         {
            contabasi[j].resize(4);

            if ((targetfoundbit[j] & guidesbit[guidafound][j]) == 0)
            {
               if (targetfound[j] == 'A')
                  contabasi[j][0]++;
               else if (targetfound[j] == 'T')
                  contabasi[j][3]++;
               else if (targetfound[j] == 'G')
                  contabasi[j][2]++;
               else if (targetfound[j] == 'C')
                  contabasi[j][1]++;
               else if (targetfound[j] == 'R')
               {
                  contabasi[j][0] = contabasi[j][0] + 0.5;
                  contabasi[j][2] = contabasi[j][2] + 0.5;
               }
               else if (targetfound[j] == 'Y')
               {
                  contabasi[j][1] = contabasi[j][1] + 0.5;
                  contabasi[j][3] = contabasi[j][3] + 0.5;
               }
               else if (targetfound[j] == 'S')
               {
                  contabasi[j][2] = contabasi[j][2] + 0.5;
                  contabasi[j][1] = contabasi[j][1] + 0.5;
               }
               else if (targetfound[j] == 'W')
               {
                  contabasi[j][0] = contabasi[j][0] + 0.5;
                  contabasi[j][3] = contabasi[j][3] + 0.5;
               }
               else if (targetfound[j] == 'K')
               {
                  contabasi[j][2] = contabasi[j][2] + 0.5;
                  contabasi[j][3] = contabasi[j][3] + 0.5;
               }
               else if (targetfound[j] == 'M')
               {
                  contabasi[j][0] = contabasi[j][0] + 0.5;
                  contabasi[j][1] = contabasi[j][1] + 0.5;
               }
               else if (targetfound[j] == 'B')
               {
                  contabasi[j][1] = contabasi[j][1] + (1 / 3);
                  contabasi[j][2] = contabasi[j][2] + (1 / 3);
                  contabasi[j][3] = contabasi[j][3] + (1 / 3);
               }
               else if (targetfound[j] == 'D')
               {
                  contabasi[j][0] = contabasi[j][0] + (1 / 3);
                  contabasi[j][2] = contabasi[j][2] + (1 / 3);
                  contabasi[j][3] = contabasi[j][3] + (1 / 3);
               }
               else if (targetfound[j] == 'H')
               {
                  contabasi[j][0] = contabasi[j][0] + (1 / 3);
                  contabasi[j][1] = contabasi[j][1] + (1 / 3);
                  contabasi[j][3] = contabasi[j][3] + (1 / 3);
               }
               else if (targetfound[j] == 'V')
               {
                  contabasi[j][0] = contabasi[j][0] + (1 / 3);
                  contabasi[j][1] = contabasi[j][1] + (1 / 3);
                  contabasi[j][2] = contabasi[j][2] + (1 / 3);
               }

               targetfound[j] = tolower(targetfound[j]);
               countmissmatchpos++;
               contamissmatch[j] = 5;
               guideprofiling[guidafound][j]++;
            }
         }
      }

      if (countmissmatchpos > 0)
      {
         //aggiorno il profili della guida con OFF-TARGET E MM TOTALI
         guideprofiling[guidafound][guidelencorrected+1]++;
         guideprofiling[guidafound][guidelencorrected+2] += countmissmatchpos;
      }
      else
      {
         //aggiorno profilo con ON-TARGET
         guideprofiling[guidafound][guidelencorrected]++;
      }

      //aggiorno profilo con nuovi MM nella colonna con threshold
      guideprofiling[guidafound][guidelen + countmissmatchpos]++;

      //aggiorno extendend_profile con dati della guida
      for (int gg = 0; gg < guidelen - pamlimit; gg++)
      {
         if (contamissmatch[gg] == 5)
         {
            matrixprofiling[guidafound][countmissmatchpos][gg][4]++;
         }

         for (int dd = 0; dd < 4; dd++)
         {
            matrixprofiling[guidafound][countmissmatchpos][gg][dd] += contabasi[gg][dd];
         }
      }
      return targetfound;
   }
   else if (nowriting == 'r' && noprofile != 'p') //entro in questo if se il profilo è disattivato, conto solo MM, metto minuscole le lettere che hanno generato MM e ritorno la guida
   {
      countmissmatchpos = 0;

      if (reverse == 0)
      {
         for (j = 0; j < guidelencorrected; j++)
         {
            if ((genomebit[index + j] & guidesbit[guidafound][j]) == 0)
            {
               countmissmatchpos++;
               targetfound[j] = tolower(targetfound[j]);
            }
         }
      }
      else
      {
         for (j = 0; j < guidelencorrected; j++)
         {
            if ((targetfoundbit[j] & guidesbit[guidafound][j]) == 0)
            {
               countmissmatchpos++;
               targetfound[j] = tolower(targetfound[j]);
            }
         }
      }
      return targetfound;
   }
}

//execute guide search
void guide_searching()
{
   vector<int> respos, resneg, respos_private, resneg_private;         //inizializzazione variabili per store indici pam
   vector<int> guidepos, guideneg, guidepos_private, guideneg_private; //inizializzaione variabili per store indici guide

   int pampossize = pamindices.size();
   int pamnegsize = pamindicesreverse.size();

   //parallel region to find all targets on a selected chromosome
   #pragma omp parallel num_threads(threads) private(respos_private, guidepos_private, resneg_private, guideneg_private)
   {
      #pragma omp for schedule(static) private(j, guidecount, currentmissmatch) nowait
      for (i = 0; i < pampossize; i++) //ciclo sugli indici trovati precendetemente, PAM
      {
         for (guidecount = 0; guidecount < totalguides; guidecount++) //ciclo le guide, su tutti gli indici
         {
            currentmissmatch = 0; //missmatch di ogni guida

            for (j = 0; j < guidelencorrected; j++) //confronto lettera-lettera, guida su genoma
            {
               if ((genomebit[j + pamindices[i]] & guidesbit[guidecount][j]) == 0)
               {
                  currentmissmatch++;
               }
               if (currentmissmatch > missmatchthrestotal[guidecount]) //superata soglia missmatch massimi, esco dal ciclo della guida
               {
                  break;
               }
            }
            if (currentmissmatch <= missmatchthrestotal[guidecount]) //se rimango nella soglia, salvo dati, guida e pam che hanno generato target
            {
               respos_private.push_back(pamindices[i]);
               guidepos_private.push_back(guidecount);
            }
         }
      }

      #pragma omp for schedule(static) private(j, guidecount, currentmissmatch) nowait
      for (i = 0; i < pamnegsize; i++)
      {
         for (guidecount = 0; guidecount < totalguides; guidecount++)
         {
            currentmissmatch = 0;

            for (j = pamlimit; j < guidelen; j++)
            {
               if ((genomebit[j + pamindicesreverse[i]] & reverseguidesbit[guidecount][j]) == 0)
               {
                  currentmissmatch++;
               }
               if (currentmissmatch > missmatchthrestotal[guidecount])
               {
                  break;
               }
            }
            if (currentmissmatch <= missmatchthrestotal[guidecount])
            {
               resneg_private.push_back(pamindicesreverse[i]);
               guideneg_private.push_back(guidecount);
            }
         }
      }
      #pragma omp critical
      {
         //copy all threads values in a single value usable to write results
         respos.insert(respos.end(), respos_private.begin(), respos_private.end());
         guidepos.insert(guidepos.end(), guidepos_private.begin(), guidepos_private.end());

         resneg.insert(resneg.end(), resneg_private.begin(), resneg_private.end());
         guideneg.insert(guideneg.end(), guideneg_private.begin(), guideneg_private.end());
      }
   }


   //check delle posizioni con risultati
   for (i = 0; i < respos.size(); i++)
   {
      totalbuf += "X"; //bulge type (X stands for no bulge)
      totalbuf += "\t";
      totalbuf += guides[guidepos[i]]; //guida
      totalbuf += "\t";
      totalbuf += missmatching(genome.substr(respos[i], guidelen), respos[i], guidepos[i], 0); //calcolo della stringa e aggiornamento contatori profili
      totalbuf += "\t";
      totalbuf += chrnames; //nome chr
      totalbuf += "\t";
      totalbuf += to_string(respos[i]); //posizione del TARGET nel genoma
      totalbuf += "\t";
      totalbuf += "+"; //strand
      totalbuf += "\t";
      totalbuf += to_string(countmissmatchpos); //mismatch per target
      totalbuf += "\t";
      totalbuf += "0"; //num bulge, sempre uguale 0
      totalbuf += "\n";
   }

   if (nowriting != 'r')
   {
      totalbuf.clear();
   }
   else
   {
      results << totalbuf;
      totalbuf.clear();
   }

   for (i = 0; i < resneg.size(); i++)
   {
      totalbuf += "X"; //bulge type (X stands for no bulge)
      totalbuf += "\t";
      totalbuf += guides[guideneg[i]];
      totalbuf += "\t";
      string reversesequence = reversetarget(genome.substr(resneg[i], guidelen));
      totalbuf += missmatching(reversesequence, resneg[i], guideneg[i], 1);
      totalbuf += "\t";
      totalbuf += chrnames;
      totalbuf += "\t";
      totalbuf += to_string(resneg[i]);
      totalbuf += "\t";
      totalbuf += "-";
      totalbuf += "\t";
      totalbuf += to_string(countmissmatchpos);
      totalbuf += "\t";
      totalbuf += "0";
      totalbuf += "\n";
   }

   if (nowriting != 'r')
   {
      totalbuf.clear();
   }
   else
   {
      results << totalbuf;
      totalbuf.clear();
   }
}
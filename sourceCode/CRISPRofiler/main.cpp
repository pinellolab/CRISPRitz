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
#include <cstdlib>
#include "crispritz.h"

using namespace std;

// #define MAXS 10000
// #define MAXC 93
// #define MAXW 2000

DIR *d;
dirent *dir;
ofstream results, profile, extentedprofile;
int i, j, genlen, pamlen, guidelen, currentmissmatch, space, pamlimit, guidecount, totalguides, filenumber, threads, pampossize, pamnegsizeinputmissmatch, countmissmatchpos, guidelencorrected, inputmissmatch;
double startime, endtime, starttotal, endtotal, writestart, writend, totalallocazione, totalpam, totalguide, totallettura, totalpostprocessing, totaltimepam, totaltimeguide, totaltimereading, totaltimeguidesearch;
vector<int> missmatchthrestotal, pamindices, pamindicesreverse, totalprimenumbers;
vector<string> fileList, guides, reverseguides, writingresults, listPam;
string fastaname, pamname, guidename, exoname, introname, resultname, profilename, extendedprofilename, annotationame, chrnames, guide, pam, substring, reversesubstring, reverseguide, reversepam, line, genome, nullnucleotide, buf, totalbuf, subpos, subneg;
char nowriting, noprofile;
string writeprofile, writeextensiveprofile;
vector<bitset<4>> targetfoundbit;

//PAM automata dimensions
int g[MAXS][MAXC];
int f[MAXS];
bitset<MAXW> out[MAXS];

//profiler
vector<vector<int>> guideprofiling;                  //vettore contenente il profilo
vector<vector<vector<vector<int>>>> matrixprofiling; //vettore contenente il profilo extended
//genome and guides bit
vector<bitset<4>> genomebit;                //genoma in bit
vector<vector<bitset<4>>> guidesbit;        //guide in bit
vector<vector<bitset<4>>> reverseguidesbit; //reverse guide in bit

int main(int argc, char **argv)
{
   //start total execution time
   starttotal = omp_get_wtime();

   //assign argv variables to stream
   string resultwriting = "r";
   string profilewriting = "p";
   string profileplusresult = "t";
   nullnucleotide = "N";
   nowriting = 'n';
   noprofile = 'n';
   fastaname = argv[1];
   pamname = argv[2];
   guidename = argv[3];
   inputmissmatch = atoi(argv[4]);
   resultname = argv[5];
   resultname += ".targets.txt";
   profilename = argv[5];
   profilename += ".profile.xls";
   extendedprofilename = argv[5];
   extendedprofilename += ".extended_profile.xls";

   //setting number threads used
   threads = omp_get_max_threads();

   if (argc > 6) //controllo che l'utente abbia inserito il valore thread altrimenti metto default tutto il disponibile
   {
      threads = atoi(argv[6]);
      if (threads > omp_get_max_threads() || threads == 0)
      {
         threads = omp_get_max_threads();
      }
   }

   if (argc > 7 && (argv[7] == resultwriting)) //consenso a scrivere i result
   {
      nowriting = 'r';
      results.open(resultname);
      results << "#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge_Size\n";
   }
   else if (argc > 7 && (argv[7] == profilewriting)) //consenso a scrivere i profili (profile+extended_profile)
   {
      noprofile = 'p';
      profile.open(profilename);
      extentedprofile.open(extendedprofilename);
   }
   else if (argc > 7 && (argv[7] == profileplusresult)) //consenso a scrivere profili e result
   {
      nowriting = 'r';
      noprofile = 'p';
      results.open(resultname);
      profile.open(profilename);
      extentedprofile.open(extendedprofilename);
      results << "#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge_Size\n";
   }

   cout << "SEARCH START" << endl;

   //reading pam
   reading_pam(); //leggo il file pam dall'input

   //reading guides
   reading_guide(); //leggo il file guide dall'input

   //counting number guides

   //reverse complement of guides
   //reverse_guide();

   //generating prime numbers to avoid guides conflict
   //generateprimenumbers();

   //profiler setting
   profilersetting(); //inizializzo le matrici dei profili

   //generate all PAMs
   listPam = generatePam(pam.substr(pamlen - pamlimit, pamlimit)); //genero lista comprensiva di PAM

   buildMachine(); //generate the automata for PAM search

   cout << "USED THREADS " << threads << endl;

   //reading chromosomes and execute analysis
   reading_chromosomes(argv); //inizio la ricerca sui cromosomi

   cout << "ANALYSIS COMPLETE" << endl;

   //profiling guides
   profiler(); //scrivo la profilazione

   //close the results file
   if (argc > 7 && (argv[7] == resultwriting))
   {
      results.close(); //chiudo il file result se era aperto
   }
   else if (argc > 7 && (argv[7] == profilewriting)) //chiudo i file di profili
   {
      double start = omp_get_wtime();
      profile << writeprofile;
      extentedprofile << writeextensiveprofile;
      double end = omp_get_wtime();
      cout << "total time writing " << end - start << endl;

      profile.close();
      extentedprofile.close();
   }
   else if (argc > 7 && (argv[7] == profileplusresult)) //chiudo tutti i file aperti
   {
      results.close();

      double start = omp_get_wtime();
      profile << writeprofile;
      extentedprofile << writeextensiveprofile;
      double end = omp_get_wtime();
      cout << "total time writing " << end - start << endl;

      profile.close();
      extentedprofile.close();
   }
   //end total execution time
   endtotal = omp_get_wtime();

   cout << "total time reading " << totaltimereading << endl;
   cout << "total time guide parallelo " << totaltimeguidesearch << endl;
   cout << "total time guide " << totaltimeguide << endl;
   cout << "total time pam " << totaltimepam << endl;

   cout << "SEARCH DONE IN TIME: " << endtotal - starttotal << endl;

   return 0;
}

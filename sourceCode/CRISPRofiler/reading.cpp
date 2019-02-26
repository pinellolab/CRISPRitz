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
#include <parallel/algorithm>

using namespace std;

extern DIR *d;
extern dirent *dir;
extern int i, genlen, pamlen, space, pamlimit, totalguides, filenumber, inputmissmatch;
extern double totaltimereading;
extern vector<int> missmatchthrestotal;
extern vector<string> fileList, guides;
extern string fastaname, pamname, guidename, chrnames, pam, line, genome;

void reading_pam()
{
   ifstream pamfile(pamname);

   while (getline(pamfile, line))
   {
      transform(line.begin(), line.end(), line.begin(), ::toupper);
      space = line.find(" ");
      pam = line.substr(0, space);
      pamlimit = stoi(line.substr(space, line.length() - 1));
      pamlen = pam.length();
   }
}

void reading_guide()
{
   ifstream guidefile(guidename);

   while (getline(guidefile, line))
   {
      transform(line.begin(), line.end(), line.begin(), ::toupper);
      space = line.find_last_of("N");
      string guida = line.substr(0, space+1);
      if (guida.length() == pamlen)
      {
         guides.push_back(guida);
      }
      else
      {
         cerr<<"SOME GUIDE IS NOT THE SAME LENGTH AS PAM, PLEASE CHECK"<<endl;
         exit(0);
      }
      
      missmatchthrestotal.push_back(inputmissmatch);
   }

   totalguides = guides.size();
   guidesbitconversion();
}

void reading_chromosomes(char **argv)
{
   if (fastaname.find(".fa") != -1)
   {
      double start;
      start = omp_get_wtime();
      ifstream fasta(fastaname.substr(0,fastaname.size()-1));

      while (getline(fasta, line).good())
      {
         if (line[0] == '>')
         {
            chrnames = (line.substr(1));
         }
         else
         {
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            genome += line;
         }
      }

      double end;
      end = omp_get_wtime();

      genlen = genome.length();

      cout << "ANALYZING CHROMOSOME: " << chrnames << endl;

      cout << "reading file time " << end - start << endl;

      genomebitconversion();

      analysis();
   }
   else
   {
      //reading folder of chromosomes
      i = 0;
      d = opendir(argv[1]);
      if (d)
      {
         while ((dir = readdir(d)) != NULL)
         {
            if (!(strcmp(".", dir->d_name)))
            {
               continue;
            }
            if (!(strcmp("..", dir->d_name)))
            {
               continue;
            }
            i++;
            fileList.push_back(dir->d_name);
         }
         filenumber = fileList.size();
         //reading every fasta file in the folder and call analysis for every fasta file
         for (int i = 0; i < filenumber; i++)
         {
            double start;
            start = omp_get_wtime();

            fastaname = argv[1];
            fastaname += fileList[i];
            ifstream fasta(fastaname);

            while (getline(fasta, line).good())
            {
               if (line[0] == '>')
               {
                  chrnames = (line.substr(1));
               }
               else
               {
                  transform(line.begin(), line.end(), line.begin(), ::toupper);
                  genome += line;
               }
            }

            double end;
            end = omp_get_wtime();

            totaltimereading += end - start;

            genlen = genome.length();

            cout << "ANALYZING CHROMOSOME: " << chrnames << endl;

            cout << "reading file time " << end - start << endl;

            totaltimereading += end - start;
            start = 0;

            genomebitconversion();

            analysis();
         }
         closedir(d);
      }
   }
}

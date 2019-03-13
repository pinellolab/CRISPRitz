#include "include/crispritz.h"

using namespace std;

DIR *d;
dirent *dir;
int genlen, pamlen, guidelen, pamlimit, totalguides;
string pam, genome, chrnames;
vector<string> guides, listPam;
extern string fastaname, pamname, guidename;
extern bool pamdirection;

//generate all pam and automata construction
void pamGeneration()
{
   //generate all PAMs
   if (!pamdirection)
   {
      //genero lista comprensiva di PAM
      listPam = generatePam(pam.substr(pamlen - pamlimit, pamlimit));
   }
   else
   {
      listPam = generatePam(pam.substr(0, pamlimit));
   }

   //generate the automata for PAM search
   buildMachine();
}

void reading_pam()
{
   ifstream pamfile(pamname);
   string line;

   while (getline(pamfile, line))
   {
      transform(line.begin(), line.end(), line.begin(), ::toupper);
      int space = line.find(" ");
      pam = line.substr(0, space);
      pamlimit = stoi(line.substr(space, line.length() - 1));
      if (pamlimit < 0)
      {
         pamdirection = 1;
         pamlimit = pamlimit * -1;
      }
      pamlen = pam.length();
   }

   pamGeneration();
}

void reading_guide()
{
   ifstream guidefile(guidename);
   string line;

   while (getline(guidefile, line))
   {
      string guida;
      transform(line.begin(), line.end(), line.begin(), ::toupper);
      guida = line.substr(0, pamlen);

      if (guida.length() == pamlen)
      {
         guides.push_back(guida);
      }
      else
      {
         cerr << "SOME GUIDE IS NOT THE SAME LENGTH AS PAM, PLEASE CHECK" << endl;
         exit(0);
      }
   }

   guidelen = guides[0].size();
   totalguides = guides.size();
   guidesbitconversion();
}

void reading_chromosomes(char **argv)
{
   string line;
   vector<string> fileList;

   if (fastaname.find(".fa") != -1)
   {
      ifstream fasta(fastaname.substr(0, fastaname.size() - 1));

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

      genlen = genome.length();

      cout << "ANALYZING CHROMOSOME: " << chrnames << endl;

      genomebitconversion();

      analysis();
   }
   else
   {
      //reading folder of chromosomes
      int i = 0;
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
         int filenumber = fileList.size();
         //reading every fasta file in the folder and call analysis for every fasta file
         for (int i = 0; i < filenumber; i++)
         {
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

            genlen = genome.length();

            cout << "ANALYZING CHROMOSOME " << chrnames << " (Total progress: " << fixed << std::setprecision(1) << (100.0 * (i + 1) / filenumber) << "%)" << endl;

            genomebitconversion();

            analysis();
         }
         closedir(d);
      }
   }
}

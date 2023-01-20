#include "include/crispritz.h"

using namespace std;

DIR *d;
dirent *dir;
int genlen, pamlen, guidelen, pamlimit, totalguides;
string pam, genome, chrnames;
vector<string> guides, listPam;
extern string fastaname, pamname, guidename;
extern bool pam_at_start;

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

      // if pamlimit is negative in the file, pam is at start (3')
      if (pamlimit < 0)
      {
         pam_at_start = true;
         pamlimit = pamlimit * -1;
      }
      pamlen = pamlimit;
   }
   // extract PAM in correct position if pamlimit is negative (at start) and reverse-complement it
   if (pam_at_start)
   {
      pam = reversetarget(pam.substr(0, pamlimit));
   }
   else
   {
      pam = pam.substr(pam.length() - pamlimit);
   }
   string pam_position = "5'";
   if (pam_at_start)
   {
      pam_position = "3'";
   }
   cout << "PAM sequence is " << pam << " at " << pam_position << endl;
}

void reading_guide()
{
   ifstream guidefile(guidename);
   string line;
   guidelen = 0; // guidelen includes PAM in bps

   while (getline(guidefile, line))
   {
      transform(line.begin(), line.end(), line.begin(), ::toupper);
      // Save guide into vector of guides (guide includes PAM with Ns)
      if (pam_at_start)
      {
         guides.push_back(reversetarget(line)); // reverse-complement the guide if PAM is at start (3')
      }
      else
      {
         guides.push_back(line);
      }
      // set guidelen to max length possibile of any guide
      if (line.length() > guidelen)
         guidelen = line.length();
   }
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
      // reading folder of chromosomes
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
         // reading every fasta file in the folder and call analysis for every fasta file
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

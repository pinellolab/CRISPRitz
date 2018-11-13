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

extern DIR *d;
extern dirent *dir;
extern ofstream results;
extern int i,j,genlen,pamlen,guidelen,currentmissmatch,space,pamlimit,guidecount,totalguides,filenumber,threads,pampossize,pamnegsize,jk,inizio,fine,inputmissmatch;
extern double startime,endtime,starttotal,endtotal,writestart,writend,totalallocazione,totalpam,totalguide,totallettura,totalpostprocessing,totaltimepam,totaltimeguide,totaltimereading;
extern vector<int> missmatchthrestotal,pamindices,pamindicesreverse,totalprimenumbers;
extern vector<string> fileList,guides,reverseguides,writingresults,listPam;
extern string indel,fastaname,pamname,guidename,exoname,introname,resultname,profilename,chrnames,guide,pam,substring,reversesubstring,reverseguide,reversepam,line,genome,nullnucleotide,buf,totalbuf,subpos,subneg;
int filenumber_special;
vector<string> fileList_special;
DIR *d_special;
dirent *dir_special;

//inizializzazione pos
extern int* pamindicesgpupos;
extern int* respos;
extern int* currentmissmatchpampos;
extern int* guidefoundpos;

//inizializzazione neg
extern int* pamindicesgpuneg;
extern int* resneg;
extern int* currentmissmatchpamneg;
extern int* guidefoundneg;

//inizializzazione generale
extern char* missmatchthrespostotalgpu;
extern char* genomeintero;
extern char* guideintere;
extern char* guideinterereverse;
extern int* primenumbers;

void reading_pam()
{
    ifstream pamfile(pamname);

    while(getline(pamfile,line))
    {
        transform(line.begin(), line.end(),line.begin(), ::toupper);
        space=line.find(" ");
        pam=line.substr(0,space);
        reversepam=pam;
        pamlimit=stoi(line.substr(space,line.length()-1));
        pamlen=pam.length();
    }
}

void reading_guide()
{
    ifstream guidefile(guidename);

    while(getline(guidefile,line))
    {
        transform(line.begin(), line.end(),line.begin(), ::toupper);
        space=line.find("\r");
        string guida = line.substr(0,space);
        if(guida.length()==pamlen)
        {
            guides.push_back(guida);
        }
        //guides.push_back(line);
        //reverseguides.push_back(line.substr(0,space));
        //missmatchthrestotal.push_back(stoi(line.substr(space,line.length()-1)));   
        missmatchthrestotal.push_back(inputmissmatch);
    }

    totalguides=guides.size();

    guidesbitconversion();
}

void reading_chromosomes(char **argv)
{
    if(fastaname.find(".fa")!=-1)
    {
        double start;
        start=omp_get_wtime();
        ifstream fasta(fastaname);
                                
        while(getline(fasta,line).good())
        {
            if(line[0] == '>' )
            {
                chrnames=(line.substr(1));
            }
            else
            {
                transform(line.begin(), line.end(),line.begin(), ::toupper);
                genome+=line;
            }
        }
        
        double end;
        end=omp_get_wtime();
        
        genlen=genome.length();

        cout << "ANALYZING CHROMOSOME: " << chrnames << endl;

        cout<<"reading file time "<<end-start<<endl;

        genomebitconversion();

        analysis();
    }
    else
    {
        //reading folder of chromosomes
        i=0;
        d = opendir(argv[1]);
        if(d)
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
            for(int i=0;i<filenumber;i++) 
            {
                double start;
                start=omp_get_wtime();

                fastaname=argv[1];
                fastaname+=fileList[i];
                ifstream fasta(fastaname);
                                        
                while(getline(fasta,line).good())
                {
                    if(line[0] == '>' )
                    {
                        chrnames=(line.substr(1));
                    }
                    else
                    {
                        transform(line.begin(), line.end(),line.begin(), ::toupper);
                        genome+=line;

                    }
                }
                    
                double end;
                end=omp_get_wtime();

                totaltimereading+=end-start;

                genlen=genome.length();
    
                cout << "ANALYZING CHROMOSOME: " << chrnames << endl;

                cout<<"reading file time "<<end-start<<endl;

                totaltimereading+=end-start;
                start=0;

                genomebitconversion();

                analysis();
            }
        closedir(d);
        }
    }
}

// void read_special(string path)
// {
//     //reading folder of chromosomes
//     int ki=0;
//     const char *c = path.c_str();
//     d_special = opendir(c);
//     if(d_special)
//     {
//         while ((dir_special = readdir(d_special)) != NULL)
//         {
//             if (!(strcmp(".", dir_special->d_name)))
//             {
//                 continue;
//             }
//             if (!(strcmp("..", dir_special->d_name)))
//             {
//                 continue;
//             }
//             ki++;
//             fileList_special.push_back(dir_special->d_name);
//         }
//         filenumber_special = fileList_special.size();
//         //reading every fasta file in the folder and call analysis for every fasta file
//         for(int ji=0;ji<filenumber_special;ji++) 
//         {
//             double start;
//             start=omp_get_wtime();

//             if(fileList_special[ji].find(chrnames+".enriched")!=-1)
//             {

//                 path+=fileList_special[ji];
//                 ifstream fasta_special(path);
//                 genome.clear();
                                        
//                 while(getline(fasta_special,line).good())
//                 {
//                     if(line[0] == '>' )
//                     {
//                         //chrnames=(line.substr(1));
//                     }
//                     else
//                     {
//                         transform(line.begin(), line.end(),line.begin(), ::toupper);
//                         genome+=line;
//                     }
//                 }
                
//             double end;
//             end=omp_get_wtime();

//             totaltimereading+=end-start;

//             genlen=genome.length();

//             //cout << "ANALYZING CHROMOSOME: " << chrnames << endl;

//             //cout<<"reading file time "<<end-start<<endl;

//             totaltimereading+=end-start;
//             start=0;

//             genomebitconversion();

//             //analysis();

//             break;
//             }
//         }
//     closedir(d_special);
//     }
// }

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
extern int i,j,genlen,pamlen,guidelen,currentmissmatch,space,pamlimit,guidecount,totalguides,filenumber,threads,pampossize,pamnegsize,jk,inizio,fine;;
extern double totaltimeguidesearch,startime,endtime,starttotal,endtotal,writestart,writend,totalallocazione,totalpam,totalguide,totallettura,totalpostprocessing,totaltimepam,totaltimeguide,totaltimereading;
extern vector<int> missmatchthrestotal,pamindices,pamindicesreverse,totalprimenumbers;
extern vector<string> fileList,guides,reverseguides,writingresults,listPam;
extern string fastaname,pamname,guidename,chrnames,guide,pam,substring,reversesubstring,reverseguide,reversepam,line,genome,nullnucleotide,buf,totalbuf,subpos,subneg; 

//inizializzazione pos
extern int* respos;
extern int* currentmissmatchpampos;
extern int* guidefoundpos;

//inizializzazione neg
extern int* resneg;
extern int* currentmissmatchpamneg;
extern int* guidefoundneg;

//inizializzazione generale
extern unsigned long long int* primenumbers;

extern vector<vector<int>>  guideprofiling;
extern vector<vector<vector<vector<int>>>> matrixprofiling;

string writeprofile,writeextensiveprofile;


void profiler()
{
    writeprofile="GUIDE\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\t\tONT\tOFFT\tOFFT/MM\t\t0MM\t1MM\t2MM\t3MM\t4MM\t5MM\t6MM\n";

    for(i=0;i<totalguides;i++)
    {
        writeprofile+=guides[i];
        writeprofile+="\t";

        for(j=0;j<(guidelen-pamlimit);j++)
        {
            writeprofile+=to_string(guideprofiling[i][j]);
            writeprofile+="\t";
        }

        writeprofile+="\t";
        writeprofile+=to_string(guideprofiling[i][20]);

        writeprofile+="\t";
        writeprofile+=to_string(guideprofiling[i][21]);
        writeprofile+="\t";

        double offtarget=guideprofiling[i][21]; // offtarget
        double missmatch=guideprofiling[i][22]; // offtarget/missmatch
        double meanmissmatch=missmatch/offtarget;

        writeprofile+=to_string(meanmissmatch);
        writeprofile+="\t";
        writeprofile+="\t";

        for(j=23;j<(guidelen-pamlimit+10);j++)
        {
            writeprofile+=to_string(guideprofiling[i][j]);
            writeprofile+="\t";
        }

        writeprofile+="\n";

        writeextensiveprofile+=">"+guides[i];
        writeextensiveprofile+="\n";
        writeextensiveprofile+="\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tBP\tTARGETS\n";

        for(int hh=0;hh<7;hh++)
        {
            writeextensiveprofile+=to_string(hh) + " MISSMATCH";

            for(int dd=0;dd<(guidelen-pamlimit);dd++)
            {                
                writeextensiveprofile+="\t";
                writeextensiveprofile+=to_string(matrixprofiling[i][hh][dd][4]);
            }

            writeextensiveprofile+="\t";
            writeextensiveprofile+=to_string(guideprofiling[i][23+hh]);

            writeextensiveprofile+="\n";

            writeextensiveprofile+="NUCLEOTIDE A";

            for(int ff=0;ff<(guidelen-pamlimit);ff++)
            {
                writeextensiveprofile+="\t";
                writeextensiveprofile+=to_string(matrixprofiling[i][hh][ff][0]);
            }

            writeextensiveprofile+="\n";
            writeextensiveprofile+="NUCLEOTIDE C";

            for(int ff=0;ff<(guidelen-pamlimit);ff++)
            {
                writeextensiveprofile+="\t";
                writeextensiveprofile+=to_string(matrixprofiling[i][hh][ff][1]);
            }

            writeextensiveprofile+="\n";
            writeextensiveprofile+="NUCLEOTIDE G";

            for(int ff=0;ff<(guidelen-pamlimit);ff++)
            {
                writeextensiveprofile+="\t";
                writeextensiveprofile+=to_string(matrixprofiling[i][hh][ff][2]);
            }

            writeextensiveprofile+="\n";
            writeextensiveprofile+="NUCLEOTIDE T";

            for(int ff=0;ff<(guidelen-pamlimit);ff++)
            {
                writeextensiveprofile+="\t";
                writeextensiveprofile+=to_string(matrixprofiling[i][hh][ff][3]);
            }

            writeextensiveprofile+="\n";
            writeextensiveprofile+="\n";

        }
    }
}
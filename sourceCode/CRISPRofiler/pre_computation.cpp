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
extern int i,j,genlen,pamlen,guidelen,currentmissmatch,space,pamlimit,guidecount,totalguides,filenumber,threads,pampossize,pamnegsize,jk,inizio,fine;
extern double startime,endtime,starttotal,endtotal,writestart,writend,totalallocazione,totalpam,totalguide,totallettura,totalpostprocessing,totaltimepam,totaltimeguide,totaltimereading;
extern vector<int> missmatchthrestotal,pamindices,pamindicesreverse,totalprimenumbers;
extern vector<string> fileList,guides,reverseguides,writingresults,listPam;
extern string fastaname,pamname,guidename,chrnames,guide,pam,substring,reversesubstring,reverseguide,reversepam,line,genome,nullnucleotide,buf,totalbuf,subpos,subneg; 

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
extern unsigned long long int* primenumbers;
extern vector<int> guidelegenda;


//profiler
vector<vector<int>>  guideprofiling;
vector<vector<vector<vector<int>>>> matrixprofiling;

//genome and guides bit
vector<bitset<4>> genomebit;
vector<vector<bitset<4>>> guidesbit;
vector<vector<bitset<4>>> reverseguidesbit;
int guidelencorrected;

void genomebitconversion()
{
    genomebit.resize(genlen);

    #pragma omp parallel for num_threads(threads) schedule(static)
    for(i=0;i<genlen;i++)
    {
        if(genome[i]=='A')
        {
            genomebit[i]=bitset<4>(string("0001"));;
        }
        else if(genome[i]=='C')
        {
            genomebit[i]=bitset<4>(string("0010"));;
        }
        else if(genome[i]=='G')
        {
            genomebit[i]=bitset<4>(string("0100"));;
        }
        else if(genome[i]=='T')
        {
            genomebit[i]=bitset<4>(string("1000"));;
        }
        else if(genome[i]=='N')
        {
            genomebit[i]=bitset<4>(string("0000"));;
        }
        else if(genome[i]=='R')
        {
            genomebit[i]=bitset<4>(string("0101"));;
        }
        else if(genome[i]=='Y')
        {
            genomebit[i]=bitset<4>(string("1010"));;
        }
        else if(genome[i]=='S')
        {
            genomebit[i]=bitset<4>(string("0110"));;
        }
        else if(genome[i]=='W')
        {
            genomebit[i]=bitset<4>(string("1001"));;
        }
        else if(genome[i]=='K')
        {
            genomebit[i]=bitset<4>(string("1100"));;
        }
        else if(genome[i]=='M')
        {
            genomebit[i]=bitset<4>(string("0011"));;
        }
        else if(genome[i]=='B')
        {
            genomebit[i]=bitset<4>(string("1110"));;
        }
        else if(genome[i]=='D')
        {
            genomebit[i]=bitset<4>(string("1101"));;
        }
        else if(genome[i]=='H')
        {
            genomebit[i]=bitset<4>(string("1011"));;
        }
        else if(genome[i]=='V')
        {
            genomebit[i]=bitset<4>(string("0111"));;
        }
    }
}

void guidesbitconversion()
{
    guidesbit.resize(totalguides);
    reverseguidesbit.resize(totalguides);

    for(i=0;i<totalguides;i++)
    {
        guidelen=guides[i].size();
        guidelencorrected=guidelen-pamlimit;

        guidesbit[i].resize(guidelen);
        reverseguidesbit[i].resize(guidelen);
        
        for(j=0;j<guidelencorrected;j++)
        {
            if(guides[i][j]=='A')
            {
                guidesbit[i][j]=bitset<4>(string("0001"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("1000"));
            }
            else if(guides[i][j]=='C')
            {
                guidesbit[i][j]=bitset<4>(string("0010"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("0100"));
            }
            else if(guides[i][j]=='G')
            {
                
                guidesbit[i][j]=bitset<4>(string("0100"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("0010"));
            }
            else if(guides[i][j]=='T')
            {
                guidesbit[i][j]=bitset<4>(string("1000"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("0001"));
            }
            else if(guides[i][j]=='N')
            {
                guidesbit[i][j]=bitset<4>(string("0000"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("0000"));
            }
            else if(guides[i][j]=='R')
            {
                guidesbit[i][j]=bitset<4>(string("0101"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("1010"));
            }
            else if(guides[i][j]=='Y')
            {
                guidesbit[i][j]=bitset<4>(string("1010"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("0101"));
            }
            else if(guides[i][j]=='S')
            {
                guidesbit[i][j]=bitset<4>(string("0110"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("1001"));
            }
            else if(guides[i][j]=='W')
            {
                guidesbit[i][j]=bitset<4>(string("1001"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("0110"));
            }
            else if(guides[i][j]=='K')
            {
                guidesbit[i][j]=bitset<4>(string("1100"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("0011"));
            }
            else if(guides[i][j]=='M')
            {
                guidesbit[i][j]=bitset<4>(string("0011"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("1100"));
            }
            else if(guides[i][j]=='B')
            {
                guidesbit[i][j]=bitset<4>(string("1110"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("0111"));
            }
            else if(guides[i][j]=='D')
            {
                guidesbit[i][j]=bitset<4>(string("1101"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("1011"));
            }
            else if(guides[i][j]=='H')
            {
                guidesbit[i][j]=bitset<4>(string("1011"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("1101"));
            }
            else if(guides[i][j]=='V')
            {
                guidesbit[i][j]=bitset<4>(string("0111"));
                reverseguidesbit[i][guidelen-1-j]=bitset<4>(string("1110"));
            }
        }
    }
}

void profilersetting()
{
    guideprofiling.resize(totalguides);
    matrixprofiling.resize(totalguides);

    for(i=0;i<totalguides;i++)
    {
        guideprofiling[i].resize(guidelen-pamlimit+10);
        matrixprofiling[i].resize(7);

        for(int kk = 0; kk < 7; kk++)
        {
            matrixprofiling[i][kk].resize(guidelen-pamlimit);

            for(int jj=0;jj<guidelen-pamlimit;jj++)
            {
                matrixprofiling[i][kk][jj].resize(5);
            }
        }
    }
}

void resetter(int inizio)
{
    guidelegenda.clear();
    guidelegenda.resize(100);

    for(i=0;i<100;i++)
    {
        guidelegenda[i]=inizio;
        inizio++;
    }
}

void generateprimenumbers() 
{
    const int N = 100;    // number of primes
    primenumbers=new unsigned long long int[N];

    int n = 0;
    bool isPrime;

    for(i=2; n<N; i++)
    {
    // Check if i is a prime, i.e.
    // if it is divisible by one of the found primes
        isPrime = true;

        for(j=0; j<n; j++)
        {
            if(i % primenumbers[j] == 0)
            {
                isPrime = false;
                break;
            }
        }

        if(isPrime)
        {
            primenumbers[n]=i; // add to the list of primes
            n++;        // increment the number of found primes
        }
    }
}
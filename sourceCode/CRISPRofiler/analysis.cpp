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
extern double startime,endtime,starttotal,endtotal,writestart,writend,totalallocazione,totalpam,totalguide,totallettura,totalpostprocessing,totaltimepam,totaltimeguide,totaltimereading;
extern vector<int> missmatchthrestotal,pamindices,pamindicesreverse,totalprimenumbers;
extern vector<string> fileList,guides,reverseguides,writingresults,listPam;
extern string fastaname,pamname,guidename,exoname,introname,resultname,profilename,chrnames,guide,pam,substring,reversesubstring,reverseguide,reversepam,line,genome,nullnucleotide,buf,totalbuf,subpos,subneg;

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

extern vector<bitset<4>> genomebit;

void analysis()
{
    //variable reset
    pamindices.clear();
    pamindicesreverse.clear();
    buf.clear();
    pamindices.resize(genlen);
    pamindicesreverse.resize(genlen);		
    
    //setting pam and creating array
    double start;

    start=omp_get_wtime();

    listPam = generatePam(pam.substr(pamlen-pamlimit,pamlimit)); 
    string arrr[listPam.size()];  
    
    for(i=0;i<listPam.size();i++)
    {
        arrr[i]=listPam[i];
    }

    //start pam searching
    searchPam(arrr,listPam.size());

    double end;
    
    end=omp_get_wtime();

    cout<<"pam search time "<<end-start<<endl;
    totaltimepam+=end-start;

    start=omp_get_wtime();
    //start guides searching, executed number of guides times
    guide_searching();

    end=omp_get_wtime();
    
    cout<<"guide search time "<<end-start<<endl;
    totaltimeguide+=end-start;
    
    //clear the genome string for the next genome analysis
    genome.clear();
    genomebit.clear();
}
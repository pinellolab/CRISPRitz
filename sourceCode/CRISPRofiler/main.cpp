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

DIR *d;
dirent *dir;
ofstream results,profile,extentedprofile;
int i,j,genlen,pamlen,guidelen,currentmissmatch,space,pamlimit,guidecount,totalguides,filenumber,threads,pampossize,pamnegsize,jk,inizio,fine,inputmissmatch;
double startime,endtime,starttotal,endtotal,writestart,writend,totalallocazione,totalpam,totalguide,totallettura,totalpostprocessing,totaltimepam,totaltimeguide,totaltimereading,totaltimeguidesearch;
vector<int> missmatchthrestotal,pamindices,pamindicesreverse,totalprimenumbers;
vector<string> fileList,guides,reverseguides,writingresults,listPam;
string fastaname,pamname,guidename,exoname,introname,resultname,profilename,extendedprofilename,annotationame,chrnames,guide,pam,substring,reversesubstring,reverseguide,reversepam,line,genome,nullnucleotide,buf,totalbuf,subpos,subneg; 
char nowriting,noprofile;

//inizializzazione pos
int* pamindicesgpupos;
int* respos;
int* currentmissmatchpampos;
unsigned long long int* guidefoundpos;

//inizializzazione neg
int* pamindicesgpuneg;
int* resneg;
int* currentmissmatchpamneg;
unsigned long long int* guidefoundneg;

//inizializzazione generale
char* missmatchthrespostotalgpu;
char* genomeintero;
char* guideintere;
char* guideinterereverse;
unsigned long long int* primenumbers;

extern string writeprofile,writeextensiveprofile;

int main( int argc, char **argv )
{
    //start total execution time
    starttotal=omp_get_wtime();

    //assign argv variables to stream
    string resultwriting = "r";
    string profilewriting = "p";
    string profileplusresult = "t";
    nullnucleotide = "N";
    nowriting='n';
    noprofile='n';
    fastaname=argv[1];
    pamname=argv[2];
    guidename=argv[3];
    inputmissmatch=atoi(argv[4]);
    resultname=argv[5];
    resultname+=".targets.txt";
    profilename=argv[5];
    profilename+=".profile.xls";
    extendedprofilename=argv[5];
    extendedprofilename+=".extended_profile.xls";

    //setting number threads used
    threads=omp_get_max_threads();

    if(argc>6)
    {
        threads=atoi(argv[6]);
        if(threads>omp_get_max_threads() || threads==0)
        {
            threads=omp_get_max_threads();
        }
    }
    
    if(argc>7 && (argv[7] == resultwriting))
    {
        nowriting='r';
        results.open(resultname);
    }
    else if(argc>7 && (argv[7] == profilewriting))
    {
        noprofile='p';
        profile.open(profilename);
        extentedprofile.open(extendedprofilename);
    }
    else if(argc>7 && (argv[7] == profileplusresult))
    {
        nowriting='r';
        noprofile='p';
        results.open(resultname);
        profile.open(profilename);
        extentedprofile.open(extendedprofilename);
    }

    cout << "SEARCH START" << endl;
    
    //reading pam
    reading_pam();
                
    //reading guides
    reading_guide();

    //counting number guides
            
    //reverse complement of guides
    //reverse_guide();

    //generating prime numbers to avoid guides conflict
    generateprimenumbers();

    //profiler setting
    profilersetting();

    //generate all PAMs
    listPam = generatePam(pam.substr(pamlen-pamlimit,pamlimit));
    

    cout << "USED THREADS " << threads << endl;        

    //reading chromosomes and execute analysis
    reading_chromosomes(argv);

    cout<<"ANALYSIS COMPLETE"<<endl;

    //profiling guides
    profiler();

    //close the results file
    if(argc>7 && (argv[7] == resultwriting))
    {
        results.close();
    }   
    else if(argc>7 && (argv[7] == profilewriting))
    {
        double start=omp_get_wtime();
        profile << writeprofile;
        extentedprofile << writeextensiveprofile;
        double end=omp_get_wtime();
        cout<<"total time writing "<<end-start<<endl;

        profile.close();
        extentedprofile.close();
    }
    else if(argc>7 && (argv[7] == profileplusresult))
    {
        results.close();

        double start=omp_get_wtime();
        profile << writeprofile;
        extentedprofile << writeextensiveprofile;
        double end=omp_get_wtime();
        cout<<"total time writing "<<end-start<<endl;

        profile.close();
        extentedprofile.close();
    }
    //end total execution time
    endtotal=omp_get_wtime();

    cout<<"total time reading "<<totaltimereading<<endl;
    cout<<"total time guide parallelo "<<totaltimeguidesearch<<endl;
    cout<<"total time guide "<<totaltimeguide<<endl;
    cout<<"total time pam "<<totaltimepam<<endl;

    cout << "SEARCH DONE IN TIME: " << endtotal-starttotal<< endl;

    return 0;
}

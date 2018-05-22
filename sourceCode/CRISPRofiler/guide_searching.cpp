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

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

extern DIR *d;
extern dirent *dir;
extern ofstream results;
extern int i,j,genlen,pamlen,guidelen,currentmissmatch,space,pamlimit,guidecount,totalguides,filenumber,threads,pampossize,pamnegsize,jk,inizio,fine;;
extern double totaltimeguidesearch,startime,endtime,starttotal,endtotal,writestart,writend,totalallocazione,totalpam,totalguide,totallettura,totalpostprocessing,totaltimepam,totaltimeguide,totaltimereading;
extern vector<int> missmatchthrestotal,pamindices,pamindicesreverse,totalprimenumbers;
extern vector<string> fileList,guides,reverseguides,writingresults,listPam;
extern string fastaname,pamname,guidename,chrnames,guide,pam,substring,reversesubstring,reverseguide,reversepam,line,genome,nullnucleotide,buf,totalbuf,subpos,subneg;
extern char nowriting;

//inizializzazione pos
extern int* respos;
extern int* currentmissmatchpampos;
extern unsigned long long int* guidefoundpos;

//inizializzazione neg
extern int* resneg;
extern int* currentmissmatchpamneg;
extern unsigned long long int* guidefoundneg;

//inizializzazione generale
extern unsigned long long int* primenumbers;

int countmissmatchpos;
int countmissmatchneg;

extern vector<vector<int>>  guideprofiling;
extern vector<vector<vector<vector<int>>>> matrixprofiling;

vector<char> contamissmatch; //conta missmatch per BP
vector<vector<float>> contabasi; //check mm nucleotides in every base
vector<int> guidelegenda;

extern vector<bitset<4>> genomebit;
extern vector<vector<bitset<4>>> guidesbit;
extern vector<vector<bitset<4>>> reverseguidesbit;
extern int guidelencorrected;

string reversetarget(string targetfound)
{
    for(j=0;j<guidelen;j++)
    {
        if (targetfound[j] == 'A') targetfound[j] = 'T';
        else if (targetfound[j] == 'T') targetfound[j] = 'A';
        else if (targetfound[j] == 'G') targetfound[j] = 'C';
        else if (targetfound[j] == 'C') targetfound[j] = 'G';
    }

    reverse(targetfound.begin(),targetfound.end());

    return targetfound;
}

string missmatching(string targetfound,int index, int guidafound)
{
    contamissmatch.clear();
    contabasi.clear();
    countmissmatchpos=0;
    contamissmatch.resize(guidelen-pamlimit);
    contabasi.resize(guidelen-pamlimit);

    for(j=0;j<guidelencorrected;j++)
    {
        contabasi[j].resize(4);

        if(targetfound[j]!=guides[guidafound][j])
        {
            if (targetfound[j] == 'A') contabasi[j][0]++;
            else if (targetfound[j] == 'T') contabasi[j][3]++;
            else if (targetfound[j] == 'G') contabasi[j][2]++;
            else if (targetfound[j] == 'C') contabasi[j][1]++;
            else if(targetfound[j] =='R') 
            {
                contabasi[j][0]=contabasi[j][0]+0.5;
                contabasi[j][2]=contabasi[j][2]+0.5;
            }
            else if(targetfound[j] =='Y') 
            {
                contabasi[j][1]=contabasi[j][1]+0.5;
                contabasi[j][3]=contabasi[j][3]+0.5;
            }
            else if(targetfound[j] =='S') 
            {
                contabasi[j][2]=contabasi[j][2]+0.5;
                contabasi[j][1]=contabasi[j][1]+0.5;
            }
            else if(targetfound[j] =='W') 
            {
                contabasi[j][0]=contabasi[j][0]+0.5;
                contabasi[j][3]=contabasi[j][3]+0.5;
            }
            else if(targetfound[j] =='K') 
            {
                contabasi[j][2]=contabasi[j][2]+0.5;
                contabasi[j][3]=contabasi[j][3]+0.5;
            }
            else if(targetfound[j] =='M') 
            {
                contabasi[j][0]=contabasi[j][1]+0.5;
                contabasi[j][1]=contabasi[j][1]+0.5;
            }
            else if(targetfound[j] =='B') 
            {
                contabasi[j][1]=contabasi[j][1]+(1/3);
                contabasi[j][2]=contabasi[j][2]+(1/3);
                contabasi[j][3]=contabasi[j][3]+(1/3);
            }
            else if(targetfound[j] =='D') 
            {
                contabasi[j][0]=contabasi[j][0]+(1/3);
                contabasi[j][2]=contabasi[j][2]+(1/3);
                contabasi[j][3]=contabasi[j][3]+(1/3);
            }
            else if(targetfound[j] =='H') 
            {
                contabasi[j][0]=contabasi[j][0]+(1/3);
                contabasi[j][1]=contabasi[j][1]+(1/3);
                contabasi[j][3]=contabasi[j][3]+(1/3);
            }
            else if(targetfound[j] =='V') 
            {
                contabasi[j][0]=contabasi[j][0]+(1/3);
                contabasi[j][1]=contabasi[j][1]+(1/3);
                contabasi[j][2]=contabasi[j][2]+(1/3);
            }

            targetfound[j]=tolower(targetfound[j]);
            countmissmatchpos++;
            contamissmatch[j]=5;
            guideprofiling[guidafound][j]++;
        }
    }

    if(countmissmatchpos>0)
    {
        guideprofiling[guidafound][21]++;
        guideprofiling[guidafound][22]+=countmissmatchpos;
    }
    else
    {
        guideprofiling[guidafound][20]++;
    }

    guideprofiling[guidafound][23+countmissmatchpos]++;

    for(int gg=0;gg<guidelen-pamlimit;gg++)
    {
        if(contamissmatch[gg]==5)
        {
            matrixprofiling[guidafound][countmissmatchpos][gg][4]++;
        }

        for(int dd=0;dd<4;dd++)
        {
            matrixprofiling[guidafound][countmissmatchpos][gg][dd]+=contabasi[gg][dd];
        }
    }

    contamissmatch.clear();
    contabasi.clear();

    return targetfound;
}

void checkguide(unsigned long long int guidaprima)
{

    totalprimenumbers.clear();

    for(j=0;j<100;j++)
    {
        if(guidaprima%primenumbers[j] == 0)
        {
            totalprimenumbers.push_back(guidelegenda[j]);
        }
    }
}

//execute guide search
void guide_searching()
{

    //inizializzazione pos
    pampossize=pamindices.size();
    respos=new int[pampossize];
    guidefoundpos=new unsigned long long int[pampossize];

    //inizializzazione neg
    pamnegsize=pamindicesreverse.size();
    resneg=new int[pamnegsize];
    guidefoundneg=new unsigned long long int[pamnegsize];

    double newguide=totalguides;
    //parallel region to find all targets on a selected chromosome

    for(jk=0;jk<(ceil(newguide/100));jk++)
    {
        int inizio=jk*100;
        int fine=MIN(((jk+1)*100),totalguides);

        resetter(inizio);

        #pragma omp parallel for num_threads(threads)
        for(i=0;i<pampossize;i++)
        {
            respos[i]=-5;
            guidefoundpos[i]=1;
        }

        #pragma omp parallel for num_threads(threads)
        for(i=0;i<pamnegsize;i++)
        {
            resneg[i]=-5;
            guidefoundneg[i]=1;
        }

        double timeguidestart=omp_get_wtime();

        #pragma omp parallel for schedule(static) private(j,guidecount,currentmissmatch) num_threads(threads)
        for (i=0;i<pampossize;i++) //ciclo sugli indici trovati precendetemente, PAM
        {
            for(guidecount=inizio;guidecount<fine;guidecount++) //ciclo le guide, su tutti gli indici
            {
                currentmissmatch=0; //missmatch di ogni guida

                for(j=0;j<guidelencorrected;j++) //confronto lettera-lettera, guida su genoma
                {
                    if((genomebit[j+pamindices[i]]&guidesbit[guidecount][j]) == 0000)
                    {
                        currentmissmatch++;
                    }
                    if(currentmissmatch>missmatchthrestotal[guidecount]) //superata soglia missmatch massimi, esco dal ciclo della guida
                    {
                        break;
                    }
                }
                if(currentmissmatch<=missmatchthrestotal[guidecount]) //se rimango nella soglia, copio la stringa dal genoma
                {
                    respos[i]=pamindices[i];
                    guidefoundpos[i]*=primenumbers[(guidecount%100)];
                }
            }
        }

        #pragma omp parallel for schedule(static) private(j,guidecount,currentmissmatch) num_threads(threads)
        for (i=0;i<pamnegsize;i++)
        {
            for(guidecount=inizio;guidecount<fine;guidecount++)
            {
                currentmissmatch=0;

                for(j=pamlimit;j<guidelen;j++)
                {
                    if((genomebit[j+pamindicesreverse[i]]&reverseguidesbit[guidecount][j]) == 0000)
                    {
                        currentmissmatch++;
                    }
                    if(currentmissmatch>missmatchthrestotal[guidecount])
                    {
                        break;
                    }
                }
                if(currentmissmatch<=missmatchthrestotal[guidecount])
                {
                    resneg[i]=pamindicesreverse[i];
                    guidefoundneg[i]*=primenumbers[(guidecount%100)];
                }
            }
        }

        double timeguideend=omp_get_wtime();

        totaltimeguidesearch+=timeguideend-timeguidestart;

        for(i=0;i<pampossize;i++)
        {
            if(respos[i]>=0)
            {
                checkguide(guidefoundpos[i]);

                for(int kg=0;kg<totalprimenumbers.size();kg++)
                {
                    totalbuf+=guides[totalprimenumbers[kg]];
                    totalbuf+="\t";
                    totalbuf+=chrnames;
                    totalbuf+="\t";
                    totalbuf+=to_string(respos[i]);
                    totalbuf+="\t";
                    totalbuf+=missmatching(genome.substr(respos[i],guidelen),respos[i],totalprimenumbers[kg]);
                    totalbuf+="\t";
                    totalbuf+="+";
                    totalbuf+="\t";
                    totalbuf+=to_string(countmissmatchpos);
                    totalbuf+="\n";
                }
            }
        }

        for(i=0;i<pamnegsize;i++)
        {
            if(resneg[i]>=0)
            {
                checkguide(guidefoundneg[i]);

                for(int kg=0;kg<totalprimenumbers.size();kg++)
                {
                    totalbuf+=guides[totalprimenumbers[kg]];
                    totalbuf+="\t";
                    totalbuf+=chrnames;
                    totalbuf+="\t";
                    totalbuf+=to_string(resneg[i]);
                    totalbuf+="\t";
                    string reversesequence=reversetarget(genome.substr(resneg[i],guidelen));
                    totalbuf+=missmatching(reversesequence,resneg[i],totalprimenumbers[kg]);
                    totalbuf+="\t";
                    totalbuf+="-";
                    totalbuf+="\t";
                    totalbuf+=to_string(countmissmatchpos);
                    totalbuf+="\n";
                }
            }
        }

        if(nowriting != 'r')
        {
            totalbuf.clear();
        }
        else
        {
            results << totalbuf;
            totalbuf.clear();
        }
    }

    delete[] guidefoundpos;
    delete[] respos;

    delete[] guidefoundneg;
    delete[] resneg;
}
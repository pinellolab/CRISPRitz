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

#define MAXS 10000
#define MAXC 93
#define MAXW 2000

extern DIR *d;
extern dirent *dir;
extern ofstream results;
extern int i,j,genlen,pamlen,guidelen,currentmissmatch,space,pamlimit,guidecount,totalguides,filenumber,threads,pampossize,pamnegsize,jk,inizio,fine;;
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
extern int* primenumbers;


int g[MAXS][MAXC];
int f[MAXS];
bitset<MAXW> out[MAXS];
   
void buildMachine()
{
    memset(g,-1,sizeof g);
    memset(f,0,sizeof f);

    int state = 0,currState = 0,index = 0;
    string str;
    ///Building a trie, each new node gets the next number as node-name.
    for(i = 0; i<listPam.size(); i++)
    {
        str = listPam[i];
        currState = 0;

        for(j = 0; j<str.size(); j++)
        {
            index = str[j] - 33;
            if(g[currState][index] == -1)
            {
                g[currState][index] = ++state;
            }
            currState = g[currState][index];
        }
        out[currState].set(i);
        ///stores whether i'th indexed string of arr, ends at state 'currState' or not. Thus adding the string to output by using 1 bit, hhh very memory efficient.
    }
    ///Failure function
    queue<int>q;
    int s,fail;
    for(i = 0; i<MAXC; i++)
    {
        if(g[0][i] != -1)
        {
            f[g[0][i]] = 0; ///here, depth is 1
            q.push(g[0][i]);
        }
        else
        {
            g[0][i] = 0; ///Necessary in failure alg below, non-existing char back to state 0. To stop infinite loop at line 68.
        }
    }
    while(! q.empty())
    {
        s = q.front();
        q.pop();
        for(i= 0; i<MAXC; i++)
        {
            if(g[s][i] != -1)
            {
                q.push(g[s][i]);
                fail = f[s]; ///here is the perfect place to calculate failure of g[s][i],cuz here 'state:s' is (depth-1) state of 'state:g[s][i]'.
                while(g[fail][i] == -1)
                {
                    fail = f[fail];
                }
                fail = g[fail][i];
                f[g[s][i]] = fail;
                out[g[s][i]] |= out[fail]; ///merging output of the node & it's failure node.
                ///Read the paper of aho-corasick,published in 1975.
            }
        }
    }
}

void searchPam()
{
    int state;
    int index;

    int tid;
    int chunk=genlen/threads;
    int chunkfine;
    int pamlistlen=listPam.size();
    int mez_list=pamlistlen/2;

    buildMachine();

    #pragma omp parallel private(tid,j,state,index,i,chunkfine) num_threads(threads)
    {
        tid=omp_get_thread_num();
        chunkfine=((tid+1)*chunk)+(pamlimit-1);

        for(i = (tid*chunk); i < MIN(chunkfine,genlen); i++)
        {
            index = genome[i] - 33;
            while(g[state][index] == -1)   ///If non-existing state, use failure function to support automaton.
            {
                state = f[state];
            }

            state = g[state][index]; /// traverse the trie state/node for the text

            if(out[state]==0) /// if the state has 0 output
            {
                continue;
            }

            for(j = 0; j<pamlistlen; j++) ///For finding position of search strings.
            {
                if(out[state].test(j)) /// if j'th string is in the output of state, means a match is found.
                {
                    //cout<<arr[j]<<" IS MATCHED AT POSITION: "<<i-arr[j].size()+1<<endl;
                    if (j < mez_list) 
                    {
                        pamindices[i]=((i-(pamlen-1))>0) ? (i-(pamlen-1)) : 0;
                        if((i-(pamlen-1))==0)
                        {
                            pamindices[i]=-1;
                        }
                        break;
                    }
                    else
                    {
                        pamindicesreverse[i]=((i-(pamlimit-1))>0) ? (i-(pamlimit-1)) : 0;
                        if((i-(pamlimit-1))==0)
                        {
                            pamindicesreverse[i]=-1;
                        }
                        if(pamindicesreverse[i]>(genlen-guidelen))
                        {
                            pamindicesreverse[i]=0;
                        }
                        break;
                    }
                }
            }
        }
    }

    pamindices.erase(remove(pamindices.begin(), pamindices.end(), 0), pamindices.end());
    pamindices.shrink_to_fit();
    replace(pamindices.begin(),pamindices.end(),-1,0);
    
    pamindicesreverse.erase(remove(pamindicesreverse.begin(), pamindicesreverse.end(), 0), pamindicesreverse.end());
    pamindicesreverse.shrink_to_fit();
    replace(pamindicesreverse.begin(),pamindicesreverse.end(),-1,0);
}


    //Given a symbol it return a corresponding nucleotides
    string switchSymbol(char sym)
    {
        if (sym == 'A') return "ARWMDHV";
        else if (sym == 'C') return "CYSMBHV";
        else if (sym == 'G') return "GRSKBDV";
        else if (sym == 'T') return "TYWKBDH";
        else if (sym == 'R') return "ARWMDHVSKBG";	//{'A', 'G'};
        else if (sym == 'Y') return "CYSMBHVWKDT";	//{'C', 'T'};
        else if (sym == 'S') return "CYSMBHVKDRG";	//{'G', 'C'};
        else if (sym == 'W') return "ARWMDHVYKBT";	//{'A', 'T'};
        else if (sym == 'K') return "GRSKBDVYWHT";	//{'G', 'T'};
        else if (sym == 'M') return "ARWMDHVYSBC";	//{'A', 'C'};
        else if (sym == 'B') return "CYSMBHVRKDGWT";	//{'C', 'G', 'T'};
        else if (sym == 'D') return "ARWMDHVSKBGYT";	//{'A', 'G', 'T'};
        else if (sym == 'H') return "ARWMDHVYSBCKT";	//{'A', 'C', 'T'};
        else if (sym == 'V') return "ARWMDHVYSBCKG";	//{'A', 'C', 'G'};
        else if (sym == 'N') return "ACGTRYSWKMBDHV";	//{'A', 'C', 'G', 'T'};
        string str(1, sym);
        return str;
    }

    // string switchSymbol(char sym)
    // {
    //     if (sym == 'R') return "AG";	//{'A', 'G'};
    //     else if (sym == 'Y') return "CT";	//{'C', 'T'};
    //     else if (sym == 'S') return "GC";	//{'G', 'C'};
    //     else if (sym == 'W') return "AT";	//{'A', 'T'};
    //     else if (sym == 'K') return "GT";	//{'G', 'T'};
    //     else if (sym == 'M') return "AC";	//{'A', 'C'};
    //     else if (sym == 'B') return "CGT";	//{'C', 'G', 'T'};
    //     else if (sym == 'D') return "AGT";	//{'A', 'G', 'T'};
    //     else if (sym == 'H') return "ACT";	//{'A', 'C', 'T'};
    //     else if (sym == 'V') return "ACG";	//{'A', 'C', 'G'};
    //     else if (sym == 'N') return "ACGT";	//{'A', 'C', 'G', 'T'};
    //     string str(1, sym);
    //     return str;
    // }


    // Given a pam return its reverse
    string reverse(string pamInput) 
    {
        string ret = "";
        for (int nuc = 0; nuc < pamInput.length(); nuc++) {
            if (pamInput[nuc] == 'A') ret = 'T' + ret;
            else if (pamInput[nuc] == 'C') ret = 'G' + ret;
            else if (pamInput[nuc] == 'G') ret = 'C' + ret;
            else if (pamInput[nuc] == 'T') ret = 'A' + ret;
            else if (pamInput[nuc] == 'R') ret = 'Y' + ret;		
            else if (pamInput[nuc] == 'Y') ret = 'R' + ret;
            else if (pamInput[nuc] == 'S') ret = 'W' + ret;		
            else if (pamInput[nuc] == 'W') ret = 'S' + ret;		
            else if (pamInput[nuc] == 'M') ret = 'K' + ret;		
            else if (pamInput[nuc] == 'K') ret = 'M' + ret;
            else if (pamInput[nuc] == 'H') ret = 'D' + ret;
            else if (pamInput[nuc] == 'D') ret = 'H' + ret;
            else if (pamInput[nuc] == 'B') ret = 'V' + ret;
            else if (pamInput[nuc] == 'V') ret = 'B' + ret;
            else ret = pamInput[nuc] + ret;	
        }
        return ret;
    }
        
        
    // Given a list of strings return the product between strings
    vector<string> getProducts(string s[], int s_size) 
    {
        int ch;
        int combinations = 1;
        vector<string> res;
        for (i=0; i<s_size; i++)
            combinations *= s[i].size();	
    
        for (i=0; i<s_size; i++) {
            string cur = s[i];	
            int div = combinations / cur.length();
            int count = 0;
            for (ch=0; ch<cur.length(); ch++) {
                for (int len=0; len<div; len++) {
                    if (i==0) {
                        res.push_back(string(cur.substr(ch, 1)));
                    } else {
                        string tmp = res.at(count);
                        tmp.append(string(cur.substr(ch,1)));
                        res.at(count) = tmp;
                    }
                    count++;
                }
                if ((ch == cur.length()-1) && (count <= res.size()-1) && i>0)
                    ch = -1;
            }
            combinations = div;
        }
        return res;
    }
        
    // Given a pam and a automaton it fill the automaton with each pam possible
    vector<string> generatePam(string pamInput) 
    {
        string pamSup = pamInput;						// copy the input pam
        vector<string> pam_vector;						// vector of pam
        vector<string> outPam;							// vector of pam
        string nucleotides_list[pamSup.length()];		// list of nucleotides of the pam
        for (int v=0; v<2; v++) {
            for (j = 0; j<pamSup.length(); j++)	// switch the symbols to nucleotides
            {		
                nucleotides_list[j] = switchSymbol(pamSup[j]);
            }
            pam_vector = getProducts(nucleotides_list, pamSup.length());
            outPam.insert(outPam.end(), pam_vector.begin(), pam_vector.end());
            pamSup = reverse(pamSup);
        }
        return outPam;
    }
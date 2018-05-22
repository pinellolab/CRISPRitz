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


// Max number of states in the matching machine.
// Should be equal to the sum of the length of all keywords.
const int MAXS = 500;

// Maximum number of characters in input alphabet
const int MAXC = 26;

// OUTPUT FUNCTION IS IMPLEMENTED USING out[]
// Bit i in this mask is one if the word with index i
// appears when the machine enters this state.
int out[MAXS];

// FAILURE FUNCTION IS IMPLEMENTED USING f[]
short int f[MAXS];

// GOTO FUNCTION (OR TRIE) IS IMPLEMENTED USING g[][]
short int g[MAXS][MAXC];

// Builds the string matching machine.
// arr -   array of words. The index of each keyword is important:
//         "out[state] & (1 << i)" is > 0 if we just found word[i]
//         in the text.
// Returns the number of states that the built machine has.
// States are numbered 0 up to the return value - 1, inclusive.

    int buildMatchingMachine(string arr[], int k) 
    {
    
        // Initialize all values in output function as 0.
        memset(out, 0, sizeof out);

        // Initialize all values in goto function as -1.
        memset(g, -1, sizeof g);

        // Initially, we just have the 0 state
        int states = 1;

        // Construct values for goto function, i.e., fill g[][]
        // This is same as building a Trie for arr[]
        for (int i = 0; i < k; ++i) 
        {
            const string &word = arr[i];
            short int currentState = 0;

            // Insert all characters of current word in arr[]
            for (int j = 0; j < word.size(); ++j) 
            {
                short int ch = word[j] - 'A';

                // Allocate a new node (create a new state) if a
                // node for ch doesn't exist.
                if (g[currentState][ch] == -1)
                    g[currentState][ch] = states++;

                currentState = g[currentState][ch];
            }

            // Add current word in output function
            out[currentState] |= (1 << i);
        }

        // For all characters which don't have an edge from
        // root (or state 0) in Trie, add a goto edge to state
        // 0 itself
        for (short int ch = 0; ch < MAXC; ++ch)
        {
            if (g[0][ch] == -1)
            {
                g[0][ch] = 0;
            }
        }

        // Now, let's build the failure function

        // Initialize values in fail function
        memset(f, -1, sizeof f);

        // Failure function is computed in breadth first order
        // using a queue
        queue<int> q;

        // Iterate over every possible input
        for (short int ch = 0; ch < MAXC; ++ch) 
        {
            // All nodes of depth 1 have failure function value
            // as 0. For example, in above diagram we move to 0
            // from states 1 and 3.
            if (g[0][ch] != 0) 
            {
                f[g[0][ch]] = 0;
                q.push(g[0][ch]);
            }
        }

        // Now queue has states 1 and 3
        while (q.size())
        {
            // Remove the front state from queue
            int state = q.front();
            q.pop();

            // For the removed state, find failure function for
            // all those characters for which goto function is
            // not defined.
            
            bitset<5> y(MAXC);
            for (short int ch = 0; ch <= MAXC; ch++) 
            {
                bitset<5> x(ch);
                if (x==y) break;
                
                // If goto function is defined for character 'ch'
                // and 'state'
                if (g[state][ch] != -1)
                {
                    // Find failure state of removed state
                    short int failure = f[state];

                    // Find the deepest node labeled by proper
                    // suffix of string from root to current
                    // state.
                    while (g[failure][ch] == -1)
                            failure = f[failure];

                    failure = g[failure][ch];
                    f[g[state][ch]] = failure;

                    // Merge output values
                    out[g[state][ch]] |= out[failure];

                    // Insert the next level node (of Trie) in Queue
                    q.push(g[state][ch]);
                }
            }
        }
        return states;
    }

    //INPUT: indices, indices reverse, vector of pams, length of vector of pams, genome
    void searchPam(string arr[], int k) 
    {
        int l = genlen;
        int mez_K = k/2;
        
        int currentState = 0;		// Initialize current state
        int ch;	
    
        // Build machine with goto, failure and output functions
        buildMatchingMachine(arr, k);
    
    
        // Traverse the text through the nuilt machine to find all occurrences of words in arr[]
        #pragma omp parallel for num_threads(threads) schedule(static) private(j, currentState, ch)
        for (i = 0; i < l; i++) 
        {
            // If goto is not defined, use failure function
            ch = genome[i] - 'A';
            while (g[currentState][ch] == -1)
                currentState = f[currentState];
            currentState = g[currentState][genome[i] - 'A'];
            
            // If match not found, move to next state
            if (out[currentState] == 0)
                continue;
            
            // Match found, print all matching words of arr[] using output function.
            for (j = 0; j < k; j++) 
            {
                if (out[currentState] & (1 << j)) 
                {	
                    if (j < mez_K) 
                    {
                        pamindices[i]=((i-(pamlen-1))>0) ? (i-(pamlen-1)) : 0;
                        if((i-(pamlen-1))==0)
                        {
                            pamindices[i]=-1;
                        }
                        break; 
                    }
                    pamindicesreverse[i]=((i-(pamlimit-1))>0) ? (i-(pamlimit-1)) : 0;
                    if((i-(pamlimit-1))==0)
                    {
                        pamindicesreverse[i]=-1;
                    }
                    if(pamindicesreverse[i]>(genome.size()-guidelen))
                    {
                        pamindicesreverse[i]=0;
                    }
                    break;
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
        
    // Given a symbol it return a corresponding nucleotides
    string switchSymbol(char sym)
    {
        if (sym == 'R') return "AG";	//{'A', 'G'};
        else if (sym == 'Y') return "CT";	//{'C', 'T'};
        else if (sym == 'S') return "GC";	//{'G', 'C'};
        else if (sym == 'W') return "AT";	//{'A', 'T'};
        else if (sym == 'K') return "GT";	//{'G', 'T'};
        else if (sym == 'M') return "AC";	//{'A', 'C'};
        else if (sym == 'B') return "CGT";	//{'C', 'G', 'T'};
        else if (sym == 'D') return "AGT";	//{'A', 'G', 'T'};
        else if (sym == 'H') return "ACT";	//{'A', 'C', 'T'};
        else if (sym == 'V') return "ACG";	//{'A', 'C', 'G'};
        else if (sym == 'N') return "ACGT";	//{'A', 'C', 'G', 'T'};
        string str(1, sym);
        return str;
    }
        
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
        unsigned short int ch;
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
            for (j = 0; j<pamSup.length(); j++)			// switch the symbols to nucleotides
                nucleotides_list[j] = switchSymbol(pamSup[j]);
            pam_vector = getProducts(nucleotides_list, pamSup.length());
            outPam.insert(outPam.end(), pam_vector.begin(), pam_vector.end());
            pamSup = reverse(pamSup);
        }
        
        return outPam;
    }
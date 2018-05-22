// TAKED FROM http://www.geeksforgeeks.org/aho-corasick-algorithm-pattern-searching/

// C++ program for implementation of Aho Corasick algorithm for string matching
using namespace std;
#include <bits/stdc++.h>
#include <string>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <cmath>

#include <bitset>


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
int buildMatchingMachine(string arr[], int k) {
	
    // Initialize all values in output function as 0.
    memset(out, 0, sizeof out);
 
    // Initialize all values in goto function as -1.
    memset(g, -1, sizeof g);
 
    // Initially, we just have the 0 state
    int states = 1;
 
    // Construct values for goto function, i.e., fill g[][]
    // This is same as building a Trie for arr[]
    for (int i = 0; i < k; ++i) {
        const string &word = arr[i];
        short int currentState = 0;
 
        // Insert all characters of current word in arr[]
        for (int j = 0; j < word.size(); ++j) {
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
        if (g[0][ch] == -1)
            g[0][ch] = 0;
 
    // Now, let's build the failure function
 
    // Initialize values in fail function
    memset(f, -1, sizeof f);
 
    // Failure function is computed in breadth first order
    // using a queue
    queue<int> q;
 
     // Iterate over every possible input
    for (short int ch = 0; ch < MAXC; ++ch) {
        // All nodes of depth 1 have failure function value
        // as 0. For example, in above diagram we move to 0
        // from states 1 and 3.
        if (g[0][ch] != 0) {
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
        for (short int ch = 0; ch <= MAXC; ch++) {
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
 

// This function finds all occurrences of all array words in text.
void searchWords(vector<int> &indices, string arr[], int k, string text) {

	int l = text.size();
	int nThr = omp_get_max_threads();
	indices.resize(l);
	
	int j, i = 0;
	int mez_K = k/2;
	double start, end;
    int currentState = 0;		// Initialize current state
	int ch;
	
	// Build machine with goto, failure and output functions
    buildMatchingMachine(arr, k);
    
	// Traverse the text through the nuilt machine to find all occurrences of words in arr[]
	#pragma parallel omp for private(j, currentState, ch)
	for (i = 0; i < l; i++) {
		// If goto is not defined, use failure function
		ch = text[i] - 'A';
		while (g[currentState][ch] == -1)
			currentState = f[currentState];
		currentState = g[currentState][text[i] - 'A'];
		
		// If match not found, move to next state
		if (out[currentState] == 0)
			 continue;
		
		// Match found, print all matching words of arr[] using output function.
		for (j = 0; j < k; j++) {
			if (out[currentState] & (1 << j)) {	
				if (j < mez_K) {
					indices[i] = (i - arr[j].size() + 1) - 21;
					indices[i] = (indices[i] > 0) ? indices[i] : 0; 	// check the if start position of the guide is <0
					break;
				}
				if (i < l - 22)											// chek the end position of the rev guide 
					indices[i] = -(i - arr[j].size() + 2);
				break;
			}
		}
	}
	
	indices.erase(remove(indices.begin(), indices.end(), 0), indices.end());
	indices.shrink_to_fit();
}
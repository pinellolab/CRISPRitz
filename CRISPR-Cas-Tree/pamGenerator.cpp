/* ========================================================================== */
/*                                                                            */
/*   pamGenerator.cpp                                                            */
/*   (c) 2012 Author                                                          */
/*                                                                            */
/*   Description                                                              */
/*                                                                            */
/* ========================================================================== */


#include "aho-Corasick.cpp"
#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <fstream>
#include <ctype.h>
#include <omp.h>
#include <string.h>
#include <vector>


using namespace std;

// #define char NUCLEOTIDE[4] = {'A', 'C', 'G', 'T'}


// Given a symbol it return a corresponding nucleotides
string switchSymbol(char sym) {
    switch (sym){
		case 'A': return "ARWMDHV"; break;
		case 'C': return "CYSMBHV"; break;
		case 'G': return "GRSKBDV"; break;
		case 'T': return "TYWKBDH"; break;
		case 'R': return "ARWMDHVSKBG"; break;
		case 'Y': return "CYSMBHVWKDT"; break;
		case 'S': return "CYSMBHVKDRG"; break;
		case 'W': return "ARWMDHVYKBT"; break;
		case 'K': return "GRSKBDVYWHT"; break;
		case 'M': return "ARWMDHVYSBC"; break;
		case 'B': return "CYSMBHVRKDGWT"; break;
		case 'D': return "ARWMDHVSKBGYT"; break;
		case 'H': return "ARWMDHVYSBCKT"; break;
		case 'V': return "ARWMDHVYSBCKG"; break;
		case 'N': return "ACGTRYSWKMBDHV"; break;
		default: cerr << "The symbol is not an IUPAC nucleotide" << endl; break;
	}
	string str(1, sym);
    return str;
}

/*
// Given a symbol it return a corresponding nucleotides //OLD
string switchSymbol(char sym) {
    if (sym == 'R') return "AG";		//{'A', 'G'};
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
*/
// Given a pam return its reverse
string reversenuc(string pam) {
    string ret = "";
    for (int nuc = 0; nuc < pam.length(); nuc++) {
        if (pam[nuc] == 'A') ret = 'T' + ret;
        else if (pam[nuc] == 'C') ret = 'G' + ret;
        else if (pam[nuc] == 'G') ret = 'C' + ret;
        else if (pam[nuc] == 'T') ret = 'A' + ret;
		else if (pam[nuc] == 'R') ret = 'Y' + ret;		
		else if (pam[nuc] == 'Y') ret = 'R' + ret;		
		else if (pam[nuc] == 'M') ret = 'K' + ret;		
		else if (pam[nuc] == 'K') ret = 'M' + ret;
		else if (pam[nuc] == 'H') ret = 'D' + ret;
		else if (pam[nuc] == 'D') ret = 'H' + ret;
		else if (pam[nuc] == 'B') ret = 'V' + ret;
		else if (pam[nuc] == 'V') ret = 'B' + ret;
		else ret = pam[nuc] + ret;	
	}
	return ret;
}


		
vector<string> getProducts(string s[], int s_size) {
	int combinations = 1;
	vector<string> res;
	for (unsigned int i=0; i<s_size; i++)
		combinations *= s[i].size();	

	for (unsigned int i=0; i<s_size; i++) {
		string cur = s[i];	
		int div = combinations / cur.length();
		int count = 0;
		for (unsigned int ch=0; ch<cur.length(); ch++) {
			for (int len=0; len<div; len++) {
				if (i==0) {
					res.push_back(string(cur.substr(ch, 1)));
				} else {
					string tmp = res[count];
					tmp.append(string(cur.substr(ch,1)));
					res[count] = tmp;
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
vector<string> generatePam(string pamInput) {
	string pam = pamInput;						// copy the input pam
	vector<string> pam_vector;					// vector of pam
	vector<string> out;							// vector of pam output
	string nucleotides_list[pam.length()];		// list of nucleotides of the pam
	for (int i = 0; i<2; i++) {
        for (int j = 0; j<pam.length(); j++)							// switch the symbols to nucleotides
			nucleotides_list[j] = switchSymbol(pam[j]);
		pam_vector = getProducts(nucleotides_list, pam.length());		// produce all possible pam
		out.insert(out.end(), pam_vector.begin(), pam_vector.end());
        pam = reversenuc(pam);
	}
	return out;
}
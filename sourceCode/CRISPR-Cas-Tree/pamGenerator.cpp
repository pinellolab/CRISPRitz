/* ========================================================================== */
/*                                                                            */
/*   pamGenerator.cpp                                                            */
/*   (c) 2012 Author                                                          */
/*                                                                            */
/*   Description                                                              */
/*                                                                            */
/* ========================================================================== */

// #include "aho-Corasick.cpp"
#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <fstream>
#include <ctype.h>
#include <omp.h>
#include <string.h>
#include <vector>
#include <bitset>

using namespace std;

// #define char NUCLEOTIDE[4] = {'A', 'C', 'G', 'T'}

vector<bitset<4>> genome_bit_conversion(string genome) //converto il genoma dal fasta alla versione bit
{
	vector<bitset<4>> genomeBit;
	genomeBit.clear();
	genomeBit.resize(genome.length());

#pragma omp parallel for schedule(static)
	for (int i = 0; i < genome.length(); ++i)
	{
		if (genome[i] == 'A')
		{
			genomeBit[i] = bitset<4>(string("0001"));
		}
		else if (genome[i] == 'C')
		{
			genomeBit[i] = bitset<4>(string("0010"));
		}
		else if (genome[i] == 'G')
		{
			genomeBit[i] = bitset<4>(string("0100"));
		}
		else if (genome[i] == 'T')
		{
			genomeBit[i] = bitset<4>(string("1000"));
		}
		else if (genome[i] == 'N')
		{
			genomeBit[i] = bitset<4>(string("0000"));
		}
		else if (genome[i] == 'R')
		{
			genomeBit[i] = bitset<4>(string("0101"));
		}
		else if (genome[i] == 'Y')
		{
			genomeBit[i] = bitset<4>(string("1010"));
		}
		else if (genome[i] == 'S')
		{
			genomeBit[i] = bitset<4>(string("0110"));
		}
		else if (genome[i] == 'W')
		{
			genomeBit[i] = bitset<4>(string("1001"));
		}
		else if (genome[i] == 'K')
		{
			genomeBit[i] = bitset<4>(string("1100"));
		}
		else if (genome[i] == 'M')
		{
			genomeBit[i] = bitset<4>(string("0011"));
		}
		else if (genome[i] == 'B')
		{
			genomeBit[i] = bitset<4>(string("1110"));
		}
		else if (genome[i] == 'D')
		{
			genomeBit[i] = bitset<4>(string("1101"));
		}
		else if (genome[i] == 'H')
		{
			genomeBit[i] = bitset<4>(string("1011"));
		}
		else if (genome[i] == 'V')
		{
			genomeBit[i] = bitset<4>(string("0111"));
		}
	}
	return genomeBit;
}

vector<bitset<4>> pam_bit_conversion(string PAM) //converto la pam in input da nt alla versione bit
{
	vector<bitset<4>> pam_bit;
	pam_bit.clear();
	pam_bit.resize(PAM.length());

	// #pragma omp parallel for schedule(static)
	for (int i = 0; i < PAM.length(); ++i)
	{
		if (PAM[i] == 'A')
		{
			pam_bit[i] = bitset<4>(string("0001"));
		}
		else if (PAM[i] == 'C')
		{
			pam_bit[i] = bitset<4>(string("0010"));
		}
		else if (PAM[i] == 'G')
		{
			pam_bit[i] = bitset<4>(string("0100"));
		}
		else if (PAM[i] == 'T')
		{
			pam_bit[i] = bitset<4>(string("1000"));
		}
		else if (PAM[i] == 'N')
		{
			pam_bit[i] = bitset<4>(string("1111"));
		}
		else if (PAM[i] == 'R')
		{
			pam_bit[i] = bitset<4>(string("0101"));
		}
		else if (PAM[i] == 'Y')
		{
			pam_bit[i] = bitset<4>(string("1010"));
		}
		else if (PAM[i] == 'S')
		{
			pam_bit[i] = bitset<4>(string("0110"));
		}
		else if (PAM[i] == 'W')
		{
			pam_bit[i] = bitset<4>(string("1001"));
		}
		else if (PAM[i] == 'K')
		{
			pam_bit[i] = bitset<4>(string("1100"));
		}
		else if (PAM[i] == 'M')
		{
			pam_bit[i] = bitset<4>(string("0011"));
		}
		else if (PAM[i] == 'B')
		{
			pam_bit[i] = bitset<4>(string("1110"));
		}
		else if (PAM[i] == 'D')
		{
			pam_bit[i] = bitset<4>(string("1101"));
		}
		else if (PAM[i] == 'H')
		{
			pam_bit[i] = bitset<4>(string("1011"));
		}
		else if (PAM[i] == 'V')
		{
			pam_bit[i] = bitset<4>(string("0111"));
		}
	}
	return pam_bit;
}

// Given a symbol it return a corresponding nucleotides
// string switchSymbol(char sym)
// {
// 	switch (sym)
// 	{
// 	case 'A':
// 		return "ARWMDHV";
// 		break;
// 	case 'C':
// 		return "CYSMBHV";
// 		break;
// 	case 'G':
// 		return "GRSKBDV";
// 		break;
// 	case 'T':
// 		return "TYWKBDH";
// 		break;
// 	case 'R':
// 		return "ARWMDHVSKBG";
// 		break;
// 	case 'Y':
// 		return "CYSMBHVWKDT";
// 		break;
// 	case 'S':
// 		return "CYSMBHVKDRG";
// 		break;
// 	case 'W':
// 		return "ARWMDHVYKBT";
// 		break;
// 	case 'K':
// 		return "GRSKBDVYWHT";
// 		break;
// 	case 'M':
// 		return "ARWMDHVYSBC";
// 		break;
// 	case 'B':
// 		return "CYSMBHVRKDGWT";
// 		break;
// 	case 'D':
// 		return "ARWMDHVSKBGYT";
// 		break;
// 	case 'H':
// 		return "ARWMDHVYSBCKT";
// 		break;
// 	case 'V':
// 		return "ARWMDHVYSBCKG";
// 		break;
// 	case 'N':
// 		return "ACGTRYSWKMBDHV";
// 		break;
// 	default:
// 		cerr << "The symbol is not an IUPAC nucleotide" << endl;
// 		break;
// 	}
// 	string str(1, sym);
// 	return str;
// }

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
string reversenuc(string pam)
{
	string ret = "";
	for (int nuc = 0; nuc < pam.length(); nuc++)
	{
		if (pam[nuc] == 'A')
			ret = 'T' + ret;
		else if (pam[nuc] == 'C')
			ret = 'G' + ret;
		else if (pam[nuc] == 'G')
			ret = 'C' + ret;
		else if (pam[nuc] == 'T')
			ret = 'A' + ret;
		else if (pam[nuc] == 'R')
			ret = 'Y' + ret;
		else if (pam[nuc] == 'Y')
			ret = 'R' + ret;
		else if (pam[nuc] == 'M')
			ret = 'K' + ret;
		else if (pam[nuc] == 'K')
			ret = 'M' + ret;
		else if (pam[nuc] == 'H')
			ret = 'D' + ret;
		else if (pam[nuc] == 'D')
			ret = 'H' + ret;
		else if (pam[nuc] == 'B')
			ret = 'V' + ret;
		else if (pam[nuc] == 'V')
			ret = 'B' + ret;
		else
			ret = pam[nuc] + ret;
	}
	return ret;
}

// vector<string> getProducts(string s[], int s_size)
// {
// 	int combinations = 1;
// 	vector<string> res;
// 	for (unsigned int i = 0; i < s_size; i++)
// 		combinations *= s[i].size();

// 	for (unsigned int i = 0; i < s_size; i++)
// 	{
// 		string cur = s[i];
// 		int div = combinations / cur.length();
// 		int count = 0;
// 		for (unsigned int ch = 0; ch < cur.length(); ch++)
// 		{
// 			for (int len = 0; len < div; len++)
// 			{
// 				if (i == 0)
// 				{
// 					res.push_back(string(cur.substr(ch, 1)));
// 				}
// 				else
// 				{
// 					string tmp = res[count];
// 					tmp.append(string(cur.substr(ch, 1)));
// 					res[count] = tmp;
// 				}
// 				count++;
// 			}
// 			if ((ch == cur.length() - 1) && (count <= res.size() - 1) && i > 0)
// 				ch = -1;
// 		}
// 		combinations = div;
// 	}
// 	return res;
// }

vector<int> searchPAMonGenome(string pam_sequence, int pam_len, string genome_sequence, int pam_limit, bool pam_at_start, int max_bulges, int max_mismatches)
{
	vector<int> indices; //to save indices for TST extraction
	vector<bitset<4>> pam_bit = pam_bit_conversion(pam_sequence);
	vector<bitset<4>> pam_bit_reverse = pam_bit_conversion(reversenuc(pam_sequence));
	vector<bitset<4>> genome_bit = genome_bit_conversion(genome_sequence);

	if (!pam_at_start) //pam al 5' quindi in fondo alla sequenza
	{
		for (int nt = 0; nt < genome_sequence.length(); ++nt)
		{
			bool found_positive = true;
			bool found_negative = true;
			int positive_mismatches = max_mismatches;
			int negative_mismatches = max_mismatches;

			for (int pam_nt = 0; pam_nt < pam_limit; ++pam_nt)
			{
				if ((genome_bit[nt + pam_nt] & pam_bit[pam_nt]) == 0)
				{
					positive_mismatches--;
					if (positive_mismatches < 0)
						found_positive = false;
				}
				if ((genome_bit[nt + pam_nt] & pam_bit_reverse[pam_nt]) == 0)
				{
					negative_mismatches--;
					if (negative_mismatches < 0)
						found_negative = false;
				}
			}
			if (found_positive)
			{
				if (((nt + pam_limit - 1) - (pam_len - 1 + max_bulges)) >= 0) //save the pam position only if possible for a guide to attach that position(avoid out of bound)
				{
					indices.push_back(((nt + pam_limit - 1) - (pam_len - 1 + max_bulges)));
				}
			}
			if (found_negative)
			{
				if ((nt <= (genome_sequence.length() - (pam_len + max_bulges)))) //same as for positive pam(out of bound problem)
				{
					indices.push_back(-nt);
				}
			}
		}
	}
	else //pam al 3' quindi in cima alla sequenza
	{
		for (int nt = 0; nt < genome_sequence.length(); ++nt)
		{
			bool found_positive = true;
			bool found_negative = true;
			int positive_mismatches = max_mismatches;
			int negative_mismatches = max_mismatches;

			for (int pam_nt = 0; pam_nt < pam_limit; ++pam_nt)
			{
				if ((genome_bit[nt + pam_nt] & pam_bit[pam_nt]) == 0)
				{
					positive_mismatches--;
					if (positive_mismatches < 0)
						found_positive = false;
				}
				if ((genome_bit[nt + pam_nt] & pam_bit_reverse[pam_nt]) == 0)
				{
					negative_mismatches--;
					if (negative_mismatches < 0)
						found_negative = false;
				}
			}
			if (found_positive)
			{
				if ((nt <= (genome_sequence.length() - (pam_len + max_bulges)))) //save the pam position only if possible for a guide to attach that position(avoid out of bound)
				{
					indices.push_back(-nt);
				}
			}
			if (found_negative)
			{
				if (((nt + pam_limit - 1) - (pam_len - 1 + max_bulges)) >= 0) //save the pam position only if possible for a guide to attach that position(avoid out of bound)
				{
					indices.push_back(((nt + pam_limit - 1) - (pam_len - 1 + max_bulges)));
				}
			}
		}
	}
	return indices;
}
// Given a pam and a automaton it fill the automaton with each pam possible
// vector<string> generatePam(string pamInput)
// {
// 	string pam = pamInput;				   // copy the input pam
// 	vector<string> pam_vector;			   // vector of pam
// 	vector<string> out;					   // vector of pam output
// 	string nucleotides_list[pam.length()]; // list of nucleotides of the pam
// 	for (int i = 0; i < 2; i++)
// 	{
// 		for (int j = 0; j < pam.length(); j++) // switch the symbols to nucleotides
// 			nucleotides_list[j] = switchSymbol(pam[j]);
// 		pam_vector = getProducts(nucleotides_list, pam.length()); // produce all possible pam
// 		out.insert(out.end(), pam_vector.begin(), pam_vector.end());
// 		pam = reversenuc(pam);
// 	}
// 	return out;
// }

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

	// #pragma omp parallel for schedule(static)
	for (int i = 0; i < genome.length(); ++i)
	{
		if (genome[i] == 'A')
		{
			genomeBit.push_back(bitset<4>(string("0001")));
		}
		else if (genome[i] == 'C')
		{
			genomeBit.push_back(bitset<4>(string("0010")));
		}
		else if (genome[i] == 'G')
		{
			genomeBit.push_back(bitset<4>(string("0100")));
		}
		else if (genome[i] == 'T')
		{
			genomeBit.push_back(bitset<4>(string("1000")));
		}
		else if (genome[i] == 'N')
		{
			genomeBit.push_back(bitset<4>(string("0000")));
		}
		else if (genome[i] == 'R')
		{
			genomeBit.push_back(bitset<4>(string("0101")));
		}
		else if (genome[i] == 'Y')
		{
			genomeBit.push_back(bitset<4>(string("1010")));
		}
		else if (genome[i] == 'S')
		{
			genomeBit.push_back(bitset<4>(string("0110")));
		}
		else if (genome[i] == 'W')
		{
			genomeBit.push_back(bitset<4>(string("1001")));
		}
		else if (genome[i] == 'K')
		{
			genomeBit.push_back(bitset<4>(string("1100")));
		}
		else if (genome[i] == 'M')
		{
			genomeBit.push_back(bitset<4>(string("0011")));
		}
		else if (genome[i] == 'B')
		{
			genomeBit.push_back(bitset<4>(string("1110")));
		}
		else if (genome[i] == 'D')
		{
			genomeBit.push_back(bitset<4>(string("1101")));
		}
		else if (genome[i] == 'H')
		{
			genomeBit.push_back(bitset<4>(string("1011")));
		}
		else if (genome[i] == 'V')
		{
			genomeBit.push_back(bitset<4>(string("0111")));
		}
		else if (genome[i] == 'Q')
		{
			//all four bases: a SNP position whose ALTs, with an A/C/G/T reference base,
			//cover A, C, G and T (the enricher writes Q there). Not N, which stays the 0000 gap.
			genomeBit.push_back(bitset<4>(string("1111")));
		}
		else
		{
			//any other byte is a gap too. Until 2.9.1 nothing was pushed here, so the
			//scan read every later position one base off per such byte.
			genomeBit.push_back(bitset<4>(string("0000")));
		}
	}
	return genomeBit;
}

vector<bitset<4>> pam_bit_conversion(string PAM) //converto la pam in input da nt alla versione bit
{
	vector<bitset<4>> pam_bit;

	// #pragma omp parallel for schedule(static)
	for (int i = 0; i < PAM.length(); ++i)
	{
		if (PAM[i] == 'A')
		{
			pam_bit.push_back(bitset<4>(string("0001")));
		}
		else if (PAM[i] == 'C')
		{
			pam_bit.push_back(bitset<4>(string("0010")));
		}
		else if (PAM[i] == 'G')
		{
			pam_bit.push_back(bitset<4>(string("0100")));
		}
		else if (PAM[i] == 'T')
		{
			pam_bit.push_back(bitset<4>(string("1000")));
		}
		else if (PAM[i] == 'N')
		{
			pam_bit.push_back(bitset<4>(string("1111")));
		}
		else if (PAM[i] == 'R')
		{
			pam_bit.push_back(bitset<4>(string("0101")));
		}
		else if (PAM[i] == 'Y')
		{
			pam_bit.push_back(bitset<4>(string("1010")));
		}
		else if (PAM[i] == 'S')
		{
			pam_bit.push_back(bitset<4>(string("0110")));
		}
		else if (PAM[i] == 'W')
		{
			pam_bit.push_back(bitset<4>(string("1001")));
		}
		else if (PAM[i] == 'K')
		{
			pam_bit.push_back(bitset<4>(string("1100")));
		}
		else if (PAM[i] == 'M')
		{
			pam_bit.push_back(bitset<4>(string("0011")));
		}
		else if (PAM[i] == 'B')
		{
			pam_bit.push_back(bitset<4>(string("1110")));
		}
		else if (PAM[i] == 'D')
		{
			pam_bit.push_back(bitset<4>(string("1101")));
		}
		else if (PAM[i] == 'H')
		{
			pam_bit.push_back(bitset<4>(string("1011")));
		}
		else if (PAM[i] == 'V')
		{
			pam_bit.push_back(bitset<4>(string("0111")));
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

// One site of an index (docs/INDEX_FORMAT.md): ws is the window start, the + strand
// coordinate of the window's leftmost base, and is negative when the window's pad or
// flank runs off the contig start; left says which end holds the flank:
//   left  (PAM-left) : [flank K][guide G][pad mb]  -- 3' PAM - strand, 5' PAM + strand
//   right (PAM-right): [pad mb][guide G][flank K]  -- 3' PAM + strand, 5' PAM - strand
struct PamSite
{
	int ws;
	bool left;
};

// whether every position of the PAM is N: such a PAM states no requirement, and an
// index built with it stores a flank that --prune-pam filters at search time
bool pamUnconstrained(const string &pam_sequence)
{
	vector<bitset<4>> pam_bit = pam_bit_conversion(pam_sequence);
	for (size_t j = 0; j < pam_bit.size(); j++)
		if (!pam_bit[j].all())
			return false;
	return true;
}

// The sites of a contig, in the order the 2.9.1 builder found them (by the start nt
// of their K flank bases; at each nt the PAM read forward first, then reversed), so
// that a contig of bases, IUPAC codes and N with no boundary site yields the same
// leaves, tree and search output as in 2.9.1. A site is stored when its guide lies
// inside the contig:
//  - dedicated PAM: the K PAM bases match (a genomic N matches no PAM letter, not even
//    N) and lie inside the contig; the scan reaches the last PAM that fits, nt =
//    length() - K (until 2.9.1 it stopped one short on a 3' PAM, and on a 5' PAM
//    bounded on pam_len, losing every reverse-strand PAM within a guide of the end);
//  - all-N PAM: every window whose flank base next to the guide is on the contig; the
//    rest of the flank may run off it (B4: --prune-pam tests only the L bases of the
//    search PAM). A window whose flank is wholly off the contig could match no PAM.
// The pad may run off the contig or into an N in both cases (B3): mainParallel.cpp
// marks how much of it an alignment may consume. Until 2.9.1 the whole pad + guide +
// flank window had to fit, so a site within max_bulges nt of a contig end was lost
// even for a perfect match, and the sign of an index was its strand, so a window
// starting at 0 was read as the other shape.
vector<PamSite> searchPAMonGenome(string pam_sequence, int pam_len, string genome_sequence, int pam_limit, bool pam_at_start, int max_bulges, int max_mismatches)
{
	vector<PamSite> sites;
	const int n = (int)genome_sequence.length();
	const int K = pam_limit;
	const int G = pam_len - pam_limit;
	//a contig shorter than a guide plus one flank base holds no site on either strand
	//(an empty or header-only FASTA is common: the enricher writes one for a
	//chromosome with no carried indel)
	if (G < 1 || K < 1 || n < G + 1)
		return sites;
	const bool unconstrained = pamUnconstrained(pam_sequence);

	//a site's flank at nt: PAM-right when the flank lies right of the guide, i.e. the
	//guide is [nt - G, nt); PAM-left when it lies left, the guide [nt + K, nt + K + G)
	auto right_ok = [&](int nt) { return nt >= G && nt <= n - 1; };
	auto left_ok = [&](int nt) { return nt + K - 1 >= 0 && nt + K + G <= n; };
	auto push_right = [&](int nt) { sites.push_back(PamSite{nt - G - max_bulges, false}); };
	auto push_left = [&](int nt) { sites.push_back(PamSite{nt, true}); };

	if (unconstrained)
	{
		//every position is a PAM, read either way; the forward read is a PAM-right
		//site for a 3' PAM and a PAM-left one for a 5' PAM
		for (int nt = -(K - 1); nt <= n - 1; nt++)
		{
			if (!pam_at_start)
			{
				if (right_ok(nt))
					push_right(nt);
				if (left_ok(nt))
					push_left(nt);
			}
			else
			{
				if (left_ok(nt))
					push_left(nt);
				if (right_ok(nt))
					push_right(nt);
			}
		}
		return sites;
	}

	if (n < K)
		return sites;
	vector<bitset<4>> pam_bit = pam_bit_conversion(pam_sequence);
	vector<bitset<4>> pam_bit_reverse = pam_bit_conversion(reversenuc(pam_sequence));
	vector<bitset<4>> genome_bit = genome_bit_conversion(genome_sequence);
	for (int nt = 0; nt <= n - K; nt++)
	{
		bool found_positive = true;
		bool found_negative = true;
		int positive_mismatches = max_mismatches;
		int negative_mismatches = max_mismatches;
		bitset<4> bit_check = 0;

		for (int pam_nt = 0; pam_nt < pam_limit; pam_nt++)
		{
			bit_check = genome_bit[nt + pam_nt] & pam_bit[pam_nt];
			if (bit_check.none())
			{
				positive_mismatches--;
				if (positive_mismatches < 0)
					found_positive = false;
			}
			bit_check = genome_bit[nt + pam_nt] & pam_bit_reverse[pam_nt];
			if (bit_check.none())
			{
				negative_mismatches--;
				if (negative_mismatches < 0)
					found_negative = false;
			}
		}
		if (!pam_at_start) //3' PAM: forward = + strand, PAM right of the guide
		{
			if (found_positive && right_ok(nt))
				push_right(nt);
			if (found_negative && left_ok(nt))
				push_left(nt);
		}
		else //5' PAM: forward = + strand, PAM left of the guide
		{
			if (found_positive && left_ok(nt))
				push_left(nt);
			if (found_negative && right_ok(nt))
				push_right(nt);
		}
	}
	return sites;
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

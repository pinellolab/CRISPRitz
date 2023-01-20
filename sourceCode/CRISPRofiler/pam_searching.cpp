#include "include/crispritz.h"
extern vector<int> pamindices, pamindicesreverse;

using namespace std;

vector<bitset<4>> genome_bit_conversion(string genome) // converto il genoma dal fasta alla versione bit
{
	vector<bitset<4>> genomeBit;
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
	}
	return genomeBit;
}

vector<bitset<4>> pam_bit_conversion(string PAM) // converto la pam in input da nt alla versione bit
{
	vector<bitset<4>> pam_bit;
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

void searchPAMonGenome(string pam_sequence, int guide_len, string genome_sequence, int pam_limit, bool pam_at_start, int max_bulges, int max_mismatches)
{
	vector<bitset<4>> pam_bit = pam_bit_conversion(pam_sequence);
	vector<bitset<4>> pam_bit_reverse = pam_bit_conversion(reversenuc(pam_sequence));
	vector<bitset<4>> genome_bit = genome_bit_conversion(genome_sequence);

	for (int nt = 0; nt < genome_sequence.length() - pam_limit; ++nt)
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
			int start_position = nt - (guide_len - pam_limit + max_bulges);
			if (start_position >= 0) // save the pam position only if possible for a guide to attach that position(avoid out of bound)
			{
				pamindices.push_back(nt); // save the starting nucleotide of the PAM, this way any guide can attach to this position
			}
		}
		if (found_negative)
		{
			int end_position = nt + (guide_len - 1) + max_bulges;
			if (end_position <= genome_sequence.length()) // same as for positive pam(out of bound problem)
			{
				pamindicesreverse.push_back(nt);
			}
		}
	}
}

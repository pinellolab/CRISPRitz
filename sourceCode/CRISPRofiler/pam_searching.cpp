// #include "include/crispritz.h"

// using namespace std;

// extern int genlen, pamlen, guidelen, pamlimit, threads;
// extern vector<int> pamindices, pamindicesreverse;
// extern vector<string> listPam;
// extern string pam, genome;
// extern bool pamdirection, variant;

// //PAM automata variables
// int g[MAXS][MAXC];
// int f[MAXS];
// vector<boost::dynamic_bitset<>> out;

// void buildMachine() //constructing aho-corasick automata
// {
// 	out.resize(listPam.size() * pamlimit);

// #pragma omp parallel for
// 	for (int ff = 0; ff < listPam.size() * pamlimit; ff++)
// 	{
// 		out[ff].resize(listPam.size());
// 	}
// 	memset(g, -1, sizeof g);
// 	memset(f, 0, sizeof f);

// 	int state = 0, currState = 0, index = 0;
// 	string str;
// 	///Building a trie, each new node gets the next number as node-name.
// 	for (int i = 0; i < listPam.size(); i++)
// 	{
// 		str = listPam[i];
// 		currState = 0;

// 		for (int j = 0; j < str.size(); j++)
// 		{
// 			index = str[j] - 33;
// 			if (g[currState][index] == -1)
// 			{
// 				g[currState][index] = ++state;
// 			}
// 			currState = g[currState][index];
// 		}
// 		out[currState].set(i);
// 		///stores whether i'th indexed string of arr, ends at state 'currState' or not. Thus adding the string to output by using 1 bit.
// 	}
// 	///Failure function
// 	queue<int> q;
// 	int s, fail;
// 	for (int i = 0; i < MAXC; i++)
// 	{
// 		if (g[0][i] != -1)
// 		{
// 			f[g[0][i]] = 0; ///here, depth is 1
// 			q.push(g[0][i]);
// 		}
// 		else
// 		{
// 			g[0][i] = 0; ///Necessary in failure alg below, non-existing char back to state 0. To stop infinite loop at line 68.
// 		}
// 	}
// 	while (!q.empty())
// 	{
// 		s = q.front();
// 		q.pop();
// 		for (int i = 0; i < MAXC; i++)
// 		{
// 			if (g[s][i] != -1)
// 			{
// 				q.push(g[s][i]);
// 				fail = f[s]; ///here is the perfect place to calculate failure of g[s][i],cuz here 'state:s' is (depth-1) state of 'state:g[s][i]'.
// 				while (g[fail][i] == -1)
// 				{
// 					fail = f[fail];
// 				}
// 				fail = g[fail][i];
// 				f[g[s][i]] = fail;
// 				out[g[s][i]] |= out[fail]; ///merging output of the node & it's failure node.
// 										   ///Read the paper of aho-corasick,published in 1975.
// 			}
// 		}
// 	}
// }

// void searchPam() //funzione che cerca le PAM nel genoma
// {
// 	//variables to maintain positions on pam search
// 	int state;
// 	int index;

// 	//variables for multithread search
// 	int tid;
// 	int chunk = 1 + ((genlen - 1) / threads);
// 	int chunkfine;
// 	int pamlistlen = listPam.size();
// 	int mez_list = pamlistlen / 2;
// 	vector<int> pamindices_private;
// 	vector<int> pamindicesreverse_private;

// #pragma omp parallel private(tid, state, index, chunkfine, pamindices_private, pamindicesreverse_private) num_threads(threads)
// 	{
// 		state = 0;
// 		index = 0;
// 		tid = omp_get_thread_num();
// 		chunkfine = ((tid + 1) * chunk) + (pamlimit - 1);

// 		for (int i = (tid * chunk); i < MIN(chunkfine, genlen); i++)
// 		{
// 			index = genome[i] - 33;
// 			while (g[state][index] == -1) ///If non-existing state, use failure function to support automaton.
// 			{
// 				state = f[state];
// 			}

// 			state = g[state][index]; /// traverse the trie state/node for the text

// 			if (!(out[state].any())) /// if the state has 0 output no match found so continue to next char on the string text
// 			{
// 				continue;
// 			}

// 			int tmp = 0;
// 			if (!pamdirection)
// 			{
// 				tmp = out[state].find_first();
// 				if (tmp < mez_list) //check what pam was found, if the first half, it's a positive pam, negative otherwise
// 				{
// 					if ((i - (pamlen - 1)) >= 0) //save the pam position only if possible for a guide to attach that position(avoid out of bound)
// 					{
// 						pamindices_private.push_back(i - (pamlen - 1));
// 					}
// 					if (out[state].find_next(tmp) != out[state].npos)
// 					{
// 						if ((i - (pamlimit - 1)) <= (genlen - guidelen)) //same as for positive pam(out of bound problem)
// 						{
// 							pamindicesreverse_private.push_back(i - (pamlimit - 1));
// 						}
// 					}
// 				}
// 				else
// 				{
// 					if ((i - (pamlimit - 1)) <= (genlen - guidelen)) //same as for positive pam(out of bound problem)
// 					{
// 						pamindicesreverse_private.push_back(i - (pamlimit - 1));
// 					}
// 				}
// 			}
// 			else
// 			{
// 				tmp = out[state].find_first();
// 				if (tmp < mez_list) //check what pam was found, if the first half, it's a positive pam, negative otherwise
// 				{
// 					if ((i - (pamlimit - 1)) <= (genlen - guidelen)) //same as for positive pam(out of bound problem)
// 					{
// 						pamindicesreverse_private.push_back(i - (pamlimit - 1));
// 					}
// 					if (out[state].find_next(tmp) != out[state].npos)
// 					{
// 						if ((i - (pamlen - 1)) >= 0) //save the pam position only if possible for a guide to attach that position(avoid out of bound)
// 						{
// 							pamindices_private.push_back(i - (pamlen - 1));
// 						}
// 					}
// 				}
// 				else
// 				{
// 					if ((i - (pamlen - 1)) >= 0) //save the pam position only if possible for a guide to attach that position(avoid out of bound)
// 					{
// 						pamindices_private.push_back(i - (pamlen - 1));
// 					}
// 				}
// 			}
// 		}
// #pragma omp critical
// 		{
// 			//copy all threads values in a single value usable to search on genome
// 			pamindices.insert(pamindices.end(), pamindices_private.begin(), pamindices_private.end());
// 			pamindicesreverse.insert(pamindicesreverse.end(), pamindicesreverse_private.begin(), pamindicesreverse_private.end());
// 		}
// 	}
// }

// //Given a symbol it return a corresponding nucleotides
// string switchSymbol_var(char sym) //creo tutte le combinazioni possibili di PAM
// {
// 	if (sym == 'A')
// 		return "ARWMDHV";
// 	else if (sym == 'C')
// 		return "CYSMBHV";
// 	else if (sym == 'G')
// 		return "GRSKBDV";
// 	else if (sym == 'T')
// 		return "TYWKBDH";
// 	else if (sym == 'R')
// 		return "ARWMDHVSKBG"; //{'A', 'G'};
// 	else if (sym == 'Y')
// 		return "CYSMBHVWKDT"; //{'C', 'T'};
// 	else if (sym == 'S')
// 		return "CYSMBHVKDRG"; //{'G', 'C'};
// 	else if (sym == 'W')
// 		return "ARWMDHVYKBT"; //{'A', 'T'};
// 	else if (sym == 'K')
// 		return "GRSKBDVYWHT"; //{'G', 'T'};
// 	else if (sym == 'M')
// 		return "ARWMDHVYSBC"; //{'A', 'C'};
// 	else if (sym == 'B')
// 		return "CYSMBHVRKDGWT"; //{'C', 'G', 'T'};
// 	else if (sym == 'D')
// 		return "ARWMDHVSKBGYT"; //{'A', 'G', 'T'};
// 	else if (sym == 'H')
// 		return "ARWMDHVYSBCKT"; //{'A', 'C', 'T'};
// 	else if (sym == 'V')
// 		return "ARWMDHVYSBCKG"; //{'A', 'C', 'G'};
// 	else if (sym == 'N')
// 		return "ACGTRYSWKMBDHV"; //{'A', 'C', 'G', 'T'};
// 	string str(1, sym);
// 	return str;
// }

// string switchSymbol_novar(char sym)
// {
// 	if (sym == 'R')
// 		return "AG"; //{'A', 'G'};
// 	else if (sym == 'Y')
// 		return "CT"; //{'C', 'T'};
// 	else if (sym == 'S')
// 		return "GC"; //{'G', 'C'};
// 	else if (sym == 'W')
// 		return "AT"; //{'A', 'T'};
// 	else if (sym == 'K')
// 		return "GT"; //{'G', 'T'};
// 	else if (sym == 'M')
// 		return "AC"; //{'A', 'C'};
// 	else if (sym == 'B')
// 		return "CGT"; //{'C', 'G', 'T'};
// 	else if (sym == 'D')
// 		return "AGT"; //{'A', 'G', 'T'};
// 	else if (sym == 'H')
// 		return "ACT"; //{'A', 'C', 'T'};
// 	else if (sym == 'V')
// 		return "ACG"; //{'A', 'C', 'G'};
// 	else if (sym == 'N')
// 		return "ACGT"; //{'A', 'C', 'G', 'T'};
// 	string str(1, sym);
// 	return str;
// }

// // Given a pam return its reverse
// string reverse(string pamInput)
// {
// 	string ret = "";
// 	for (int nuc = 0; nuc < pamInput.length(); nuc++)
// 	{
// 		if (pamInput[nuc] == 'A')
// 			ret = 'T' + ret;
// 		else if (pamInput[nuc] == 'C')
// 			ret = 'G' + ret;
// 		else if (pamInput[nuc] == 'G')
// 			ret = 'C' + ret;
// 		else if (pamInput[nuc] == 'T')
// 			ret = 'A' + ret;
// 		else if (pamInput[nuc] == 'R')
// 			ret = 'Y' + ret;
// 		else if (pamInput[nuc] == 'Y')
// 			ret = 'R' + ret;
// 		else if (pamInput[nuc] == 'S')
// 			ret = 'W' + ret;
// 		else if (pamInput[nuc] == 'W')
// 			ret = 'S' + ret;
// 		else if (pamInput[nuc] == 'M')
// 			ret = 'K' + ret;
// 		else if (pamInput[nuc] == 'K')
// 			ret = 'M' + ret;
// 		else if (pamInput[nuc] == 'H')
// 			ret = 'D' + ret;
// 		else if (pamInput[nuc] == 'D')
// 			ret = 'H' + ret;
// 		else if (pamInput[nuc] == 'B')
// 			ret = 'V' + ret;
// 		else if (pamInput[nuc] == 'V')
// 			ret = 'B' + ret;
// 		else
// 			ret = pamInput[nuc] + ret;
// 	}
// 	return ret;
// }

// // Given a list of strings return the product between strings
// vector<string> getProducts(string s[], int s_size)
// {
// 	int ch;
// 	int combinations = 1;
// 	vector<string> res;
// 	for (int i = 0; i < s_size; i++)
// 		combinations *= s[i].size();

// 	for (int i = 0; i < s_size; i++)
// 	{
// 		string cur = s[i];
// 		int div = combinations / cur.length();
// 		int count = 0;
// 		for (ch = 0; ch < cur.length(); ch++)
// 		{
// 			for (int len = 0; len < div; len++)
// 			{
// 				if (i == 0)
// 				{
// 					res.push_back(string(cur.substr(ch, 1)));
// 				}
// 				else
// 				{
// 					string tmp = res.at(count);
// 					tmp.append(string(cur.substr(ch, 1)));
// 					res.at(count) = tmp;
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

// // Given a pam and a automaton it fill the automaton with each pam possible
// vector<string> generatePam(string pamInput)
// {
// 	string pamSup = pamInput;				  // copy the input pam
// 	vector<string> pam_vector;				  // vector of pam
// 	vector<string> outPam;					  // vector of pam
// 	string nucleotides_list[pamSup.length()]; // list of nucleotides of the pam
// 	for (int v = 0; v < 2; v++)
// 	{
// 		for (int j = 0; j < pamSup.length(); j++) // switch the symbols to nucleotides
// 		{
// 			if (variant)
// 			{
// 				nucleotides_list[j] = switchSymbol_var(pamSup[j]);
// 			}
// 			else
// 			{
// 				nucleotides_list[j] = switchSymbol_novar(pamSup[j]);
// 			}
// 		}
// 		pam_vector = getProducts(nucleotides_list, pamSup.length());
// 		outPam.insert(outPam.end(), pam_vector.begin(), pam_vector.end());
// 		pamSup = reverse(pamSup);
// 	}
// 	return outPam;
// }

/* ========================================================================== */
/*                                                                            */
/*   pamGenerator.cpp                                                            */
/*   (c) 2012 Author                                                          */
/*                                                                            */
/*   Description                                                              */
/*                                                                            */
/* ========================================================================== */

// #include "aho-Corasick.cpp"
#include "include/crispritz.h"
extern vector<int> pamindices, pamindicesreverse;

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

void searchPAMonGenome(string pam_sequence, int pam_len, string genome_sequence, int pam_limit, bool pam_at_start, int max_bulges, int max_mismatches)
{
	//extract the PAM from the whole sequence
	if (!pam_at_start)
	{
		pam_sequence = pam_sequence.substr(pam_len - pam_limit, pam_len);
	}
	else
	{
		pam_sequence = pam_sequence.substr(0, pam_limit);
	}

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
					pamindices.push_back(((nt + pam_limit - 1) - (pam_len - 1 + max_bulges)));
				}
			}
			if (found_negative)
			{
				if ((nt <= (genome_sequence.length() - (pam_len + max_bulges)))) //same as for positive pam(out of bound problem)
				{
					pamindicesreverse.push_back(nt);
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
					pamindicesreverse.push_back(nt);
				}
			}
			if (found_negative)
			{
				if (((nt + pam_limit - 1) - (pam_len - 1 + max_bulges)) >= 0) //save the pam position only if possible for a guide to attach that position(avoid out of bound)
				{
					pamindices.push_back(((nt + pam_limit - 1) - (pam_len - 1 + max_bulges)));
				}
			}
		}
	}
	// return indices;
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

#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <algorithm>

#include <bitset>
#include "detailedOutput.h"
#include "convert.h"
#include <stdint.h>
#include <dirent.h> //for directory reading files
#include <array>
#include <sys/stat.h>
#include <iomanip>	 // std::setprecision
#include <sys/stat.h> //for file dimension

using namespace std;

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef struct tnode *Tptr;
typedef struct tnode
{
	char splitchar;
	bitset<4> splitchar_bit;
	int lokid, eqkid, hikid = 0;
} Tnode;

typedef struct tleaf
{
	int guideIndex;
	char *guideDNA;
	int next;
	bitset<4> *guideDNA_bit;
} Tleaf;

struct file_stats_struct
{
	string file_name;
	int file_dim;
};

Tptr tree;
string pamRNA;

vector<char *> guideRNA;					 // input guide RNA
vector<string> guideRNA_s;					 //input guide RNA as strings;
vector<vector<bitset<4>>> guideRNA_bit; //input guide RNA in bitset
int len_guide;

//global vectors for saving targets results
vector<string> vecTargetOfGuide_glb, vecInGuide_glb, bulgeType_glb;
vector<char> directions_glb;
vector<int> indices_glb, mismatches_glb, bulgeSize_glb;
vector<string> chrName_glb;
bool pam_at_start;

void readPair(vector<bitset<4>> &pairNuc_bit, char (&pairNuc)[2], char &in)
{

	if ((in & 0xF0) == 0xF0)
	{
		pairNuc_bit[0] = in >> 4;
		pairNuc[0] = '_';
		return;
	}
	else
	{
		unsigned char first = in & 0xF0; // keep only 4 left-most bits
		pairNuc_bit[0] = first >> 4;		// shift them to the right to save them in a bitset<4>
		switch (first)
		{
		case 0x10:
			pairNuc[0] = 'A';
			break;
		case 0x20:
			pairNuc[0] = 'C';
			break;
		case 0x40:
			pairNuc[0] = 'G';
			break;
		case 0x80:
			pairNuc[0] = 'T';
			break;
		//case 0xF0: pairNuc[0] = 'N'; break;
		case 0x50:
			pairNuc[0] = 'R';
			break;
		case 0xA0:
			pairNuc[0] = 'Y';
			break;
		case 0x60:
			pairNuc[0] = 'S';
			break;
		case 0x90:
			pairNuc[0] = 'W';
			break;
		case 0xC0:
			pairNuc[0] = 'K';
			break;
		case 0x30:
			pairNuc[0] = 'M';
			break;
		case 0xE0:
			pairNuc[0] = 'B';
			break;
		case 0xD0:
			pairNuc[0] = 'D';
			break;
		case 0xB0:
			pairNuc[0] = 'H';
			break;
		case 0x70:
			pairNuc[0] = 'V';
			break;
		default:
			pairNuc[0] = '0';
			break;
		}
		unsigned char second = in & 0xF; // keep only 4 right-most bits
		pairNuc_bit[1] = second;
		switch (second)
		{
		case 0x1:
			pairNuc[1] = 'A';
			break;
		case 0x2:
			pairNuc[1] = 'C';
			break;
		case 0x4:
			pairNuc[1] = 'G';
			break;
		case 0x8:
			pairNuc[1] = 'T';
			break;
		//case 0xF: pairNuc[1] = 'N'; break;
		case 0x0F:
			pairNuc[1] = '_';
			break;
		case 0x05:
			pairNuc[1] = 'R';
			break;
		case 0x0A:
			pairNuc[1] = 'Y';
			break;
		case 0x06:
			pairNuc[1] = 'S';
			break;
		case 0x09:
			pairNuc[1] = 'W';
			break;
		case 0x0C:
			pairNuc[1] = 'K';
			break;
		case 0x03:
			pairNuc[1] = 'M';
			break;
		case 0x0E:
			pairNuc[1] = 'B';
			break;
		case 0x0D:
			pairNuc[1] = 'D';
			break;
		case 0x0B:
			pairNuc[1] = 'H';
			break;
		case 0x07:
			pairNuc[1] = 'V';
			break;
		default:
			pairNuc[1] = '0';
			break;
		}
	}
}

void deSerialize(vector<Tnode> &albero, ifstream &fileTree, vector<bitset<4>> &pairNuc_bit, char (&pairNuc)[2], char &in, uint8_t &flag, int &pNode)
{

	int i = pNode;
	albero[i].splitchar_bit = pairNuc_bit[flag];
	albero[i].splitchar = pairNuc[flag];
	if (flag)
	{
		fileTree.get(in);
		readPair(pairNuc_bit, pairNuc, in);
		flag = 0;
	}
	else
		flag++;
	if (pairNuc_bit[flag] != '0')
	{
		albero[i].lokid = ++pNode;
		deSerialize(albero, fileTree, pairNuc_bit, pairNuc, in, flag, pNode);
	} // go to lokid

	if (flag)
	{
		fileTree.get(in);
		readPair(pairNuc_bit, pairNuc, in);
		flag = 0;
	}
	else
		flag++;
	if (pairNuc_bit[flag] != '0')
	{
		albero[i].hikid = ++pNode;
		deSerialize(albero, fileTree, pairNuc_bit, pairNuc, in, flag, pNode);
	} // go to hikid

	if (flag)
	{
		fileTree.get(in);
		readPair(pairNuc_bit, pairNuc, in);
		flag = 0;
	}
	else
		flag++;
	if (pairNuc_bit[flag] == '_')
	{
		flag++;
		fileTree.read((char *)&(albero[i].eqkid), sizeof(int));
	}
	else
	{
		albero[i].eqkid = ++pNode;
		deSerialize(albero, fileTree, pairNuc_bit, pairNuc, in, flag, pNode);
	} // go to eqkid
}

// Load from file by using boost
Tleaf *loadTST(string path, vector<Tnode> &albero, ifstream &fileTree, int &numNodes, int &numLeaves)
{
	Tleaf *targetOnDNA;
	fileTree.open(path, ios::in | ios::binary);

	fileTree.read((char *)&numLeaves, sizeof(int));				 // read number of leaves
	targetOnDNA = (Tleaf *)malloc(numLeaves * sizeof(Tleaf)); // initialize array of targets on DNA
	char in;
	for (int i = 0; i < numLeaves; i++)
	{																						 // fill array of targets on DNA
		fileTree.read((char *)&targetOnDNA[i].guideIndex, sizeof(int)); // read index of target on DNA

		targetOnDNA[i].guideDNA_bit = new bitset<4>[pamRNA.size()];
		targetOnDNA[i].guideDNA = new char[pamRNA.size() + 1]; // initialize PAM size
		targetOnDNA[i].guideDNA[pamRNA.size()] = '\0';
		unsigned char mask;
		int k = 0;

		fileTree.get(in);

		for (int j = pamRNA.size() - 1; j > -1; j--)
		{
			if (k == 2)
			{
				fileTree.get(in);
				k = 0;
			}

			mask = in & 0xF0;
			in <<= 4;

			switch (mask)
			{
			case 0x0:
				//cout << "Zero" ;
				break;
			case 0x10:
				targetOnDNA[i].guideDNA[j] = 'A';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0001");
				//cout << "A" ;
				break;
			case 0x20:
				targetOnDNA[i].guideDNA[j] = 'C';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0010");
				//cout << "C";
				break;
			case 0x40:
				targetOnDNA[i].guideDNA[j] = 'G';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0100");
				//cout << "G" ;
				break;
			case 0x80:
				targetOnDNA[i].guideDNA[j] = 'T';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1000");
				// cout << "T" ;
				break;
			case 0x50: //R
				targetOnDNA[i].guideDNA[j] = 'R';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0101");
				// cout << "R" ;
				break;
			case 0xA0: //Y
				targetOnDNA[i].guideDNA[j] = 'Y';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1010");
				//cout << "Y";
				break;
			case 0x60: //S
				targetOnDNA[i].guideDNA[j] = 'S';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0110");
				//cout << "S" ;
				break;
			case 0x90: //W
				targetOnDNA[i].guideDNA[j] = 'W';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1001");
				//cout << "W" ;
				break;
			case 0xC0: //K
				targetOnDNA[i].guideDNA[j] = 'K';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1100");
				//cout << "K" ;
				break;
			case 0x30: //M
				targetOnDNA[i].guideDNA[j] = 'M';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0011");
				//cout << "M" ;
				break;
			case 0xE0: //B
				targetOnDNA[i].guideDNA[j] = 'B';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1110");
				//cout << "B" ;
				break;
			case 0xD0: //D
				targetOnDNA[i].guideDNA[j] = 'D';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1101");
				//cout << "D" ;
				break;
			case 0xB0: //H
				targetOnDNA[i].guideDNA[j] = 'H';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1011");
				//cout << "H";
				break;
			case 0x70: //V
				targetOnDNA[i].guideDNA[j] = 'V';
				targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0111");
				//cout << "V";
				break;
			default:
				cerr << "Error on reading the PAM, input char is not an IUPAC nucleotide: " << mask << endl;
				break;
			}
			k++;
		}

		fileTree.get(in); // read index of next PAM with same guide

		if (in == '0')
		{
			targetOnDNA[i].next = 0;
		}
		else
		{
			fileTree.read((char *)&targetOnDNA[i].next, sizeof(int));
		}
	}

	fileTree.read((char *)&numNodes, sizeof(int)); // read number of nodes

	vector<bitset<4>> pairNuc_bit(2);
	char pairNuc[2];
	uint8_t flag = 1;

	albero.resize(numNodes);
	int pNode = 0;
	fileTree.get(in);

	readPair(pairNuc_bit, pairNuc, in);

	flag = 0;
	deSerialize(albero, fileTree, pairNuc_bit, pairNuc, in, flag, pNode); // deserialize TST

	fileTree.close();
	return targetOnDNA;
}

int bulDNA, bulRNA, mm;

bool create_target;
bool create_profile;

void saveIndices(vector<bitset<4>> &inGuide_bit, vector<bitset<4>> &targetOfGuide_bit, vector<Tnode> &p, int pos_tree, int d, int bD, int bR, int bulType,
					  Tleaf *targetOnDNA, char *inGuide, char *targetOfGuide, vector<string> &vecTargetOfGuide, vector<string> &vecInGuide, vector<string> &bulgeType, vector<char> &directions, vector<int> &indices, vector<int> &mismatches, vector<int> &bulgeSize,
					  int guideI, vector<vector<vector<int>>> &profiling, vector<vector<vector<vector<vector<int>>>>> &ext_profiling, vector<vector<vector<int>>> &pd, vector<vector<vector<int>>> &pdm, vector<vector<vector<int>>> &pr, vector<vector<vector<int>>> &prm)
{
	if (p[pos_tree].lokid > 0)
		saveIndices(inGuide_bit, targetOfGuide_bit, p, p[pos_tree].lokid, d, bD, bR, bulType, targetOnDNA,
						inGuide, targetOfGuide, vecTargetOfGuide, vecInGuide, bulgeType, directions, indices, mismatches, bulgeSize,
						guideI, profiling, ext_profiling, pd, pdm, pr, prm); // go to lokid
	if (p[pos_tree].hikid > 0)
		saveIndices(inGuide_bit, targetOfGuide_bit, p, p[pos_tree].hikid, d, bD, bR, bulType, targetOnDNA,
						inGuide, targetOfGuide, vecTargetOfGuide, vecInGuide, bulgeType, directions, indices, mismatches, bulgeSize,
						guideI, profiling, ext_profiling, pd, pdm, pr, prm); // go to hikid

	if (p[pos_tree].eqkid < 0)
	{ // node is a leaf, save index and strand

		//string g = inGuide.substr(0, len_guide + bulDNA - bD);

		string g(inGuide);
		g = g.substr(0, len_guide + bulDNA - bD);
		reverse(g.begin(), g.end());

		if (pam_at_start)
		{
			g.insert(0, pamRNA);
		}
		else
		{
			g += pamRNA;
		}

		vector<bitset<4>> g_bit(len_guide + pamRNA.size() + bulDNA);

		int j = 0;

		for (j = 0; j < (len_guide + bulDNA - bD); j++)
		{
			g_bit[j] = inGuide_bit[j];
		}

		reverse(g_bit.begin(), g_bit.begin() + len_guide + bulDNA - bD);

		int pos_to_save = 0;
		for (int i = 0; i < pamRNA.size(); i++)
		{
			g_bit[g_bit.size() - pamRNA.size() - bD + i] = bitset<4>(convertCharToBitset(pamRNA[i]));
			pos_to_save = g_bit.size() - pamRNA.size() - bD + i;
		}
		pos_to_save++;
		g_bit.resize(pos_to_save);

		int index = (p[pos_tree].eqkid + 1) * -1;

		do
		{
			//string t = targetOfGuide.substr(0, len_guide + bulDNA - bD);
			string t(targetOfGuide);
			t = t.substr(0, len_guide + bulDNA - bD);
			reverse(t.begin(), t.end());

			if (pam_at_start)
			{
				string tmp(targetOnDNA[index].guideDNA);
				reverse(tmp.begin(), tmp.end());
				t.insert(0, tmp); //targetOnDNA[index].guideDNA);
										//cout << "pam target: " << targetOnDNA[index].guideDNA << endl;
			}
			else
			{
				t += targetOnDNA[index].guideDNA;
			}
			vector<bitset<4>> t_bit(len_guide + pamRNA.size() + bulDNA);

			int i = 0;
			for (i = 0; i < (len_guide + bulDNA - bD); i++)
			{
				t_bit[i] = targetOfGuide_bit[i];
			}

			reverse(t_bit.begin(), t_bit.begin() + len_guide + bulDNA - bD);

			i = len_guide + bulDNA - bD;

			for (int k = 0; k < pamRNA.size(); k++)
			{
				t_bit[i] = targetOnDNA[index].guideDNA_bit[k];
				i++;
			}

			t_bit.resize(pos_to_save);

			//Check what type of result to create
			if (create_target)
			{

				mismatches.emplace_back(mm - d); // save mismatches
				if (bulType == 0)
				{ // NO BULGE case
					bulgeType.emplace_back("X");
					bulgeSize.emplace_back(0);
				}
				else
				{ // BULGE case
					bulgeType.emplace_back(bulType < 0 ? "RNA" : "DNA");
					bulgeSize.emplace_back(bulType < 0 ? bulRNA - bR : bulDNA - bD);
				}
				if (targetOnDNA[index].guideIndex < 0)
				{ // negative strand
					indices.emplace_back(targetOnDNA[index].guideIndex * -1);
					if (pam_at_start)
						directions.emplace_back('+');
					else
					{
						directions.emplace_back('-');
					}
				}
				else
				{ // strand positive
					indices.emplace_back(targetOnDNA[index].guideIndex + 2 - (bulDNA - bD) + (bulRNA - bR));
					if (pam_at_start)
						directions.emplace_back('-');
					else
					{
						directions.emplace_back('+');
					}
				}

				//Profiling and profiling ext
				if (create_profile)
				{
					if (bulType == 0)
						detailedOutputFast(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, profiling, ext_profiling, vecInGuide, vecTargetOfGuide);
					else if (bulType > 0)
						detailedOutputFastBulgeDNA(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, bD, bulDNA, pd, pdm, vecInGuide, vecTargetOfGuide);
					else
						detailedOutputFastBulgeRNA(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, bD, bulDNA, pr, prm, vecInGuide, vecTargetOfGuide);
				}
				else
				{
					vecInGuide.emplace_back(g);		 // save guide
					vecTargetOfGuide.emplace_back(t); // save target
				}
			}
			else
			{ //no target, only profile
				if (bulType == 0)
					detailedOutputFast(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, profiling, ext_profiling, vecInGuide, vecTargetOfGuide);
				else if (bulType > 0)
					detailedOutputFastBulgeDNA(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, bD, bulDNA, pd, pdm, vecInGuide, vecTargetOfGuide);
				else
					detailedOutputFastBulgeRNA(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, bD, bulDNA, pr, prm, vecInGuide, vecTargetOfGuide);
			}

			index = (targetOnDNA[index].next + 1) * -1;

		} while (index > -1);
	}
	else if (p[pos_tree].eqkid > 0)
		saveIndices(inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, d, bD, bR, bulType, targetOnDNA,
						inGuide, targetOfGuide, vecTargetOfGuide, vecInGuide, bulgeType, directions, indices, mismatches, bulgeSize,
						guideI, profiling, ext_profiling, pd, pdm, pr, prm); //go to eqkid
}

// CONTROLLARE LA QUESTIONE BULGE CONTEMPORANEI IN DNA ED RNA
void nearsearch(char *s, int ti, int gi, vector<bitset<4>> &inGuide_bit, vector<bitset<4>> &targetOfGuide_bit, vector<Tnode> &p, int pos_tree, vector<bitset<4>> guide_bit, int pos_in_guide, int d, int bD, int bR, bool goToLoHi, int bulType,
					 Tleaf *t, char *inGuide, char *targetOfGuide, vector<string> &vTOG, vector<string> &vIG, vector<string> &bT, vector<char> &dirs, vector<int> &ind, vector<int> &mismatches, vector<int> &bulgeSize,
					 int guideI, vector<vector<vector<int>>> &profiling, vector<vector<vector<vector<vector<int>>>>> &ext_profiling, vector<vector<vector<int>>> &pd, vector<vector<vector<int>>> &pdm, vector<vector<vector<int>>> &pr, vector<vector<vector<int>>> &prm)
{ //char *inGuide, char *targetOfGuide,

	if (p[pos_tree].lokid > 0 && goToLoHi && (d > 0 || bD > 0 || bR > 0 || *s < p[pos_tree].splitchar)) // go to lokid
		nearsearch(s, ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].lokid, guide_bit, pos_in_guide, d, bD, bR, true, bulType, t,
					  inGuide, targetOfGuide, vTOG, vIG, bT, dirs, ind, mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);

	if (p[pos_tree].hikid > 0 && goToLoHi && (d > 0 || bD > 0 || bR > 0 || *s > p[pos_tree].splitchar)) // go to hikid
		nearsearch(s, ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].hikid, guide_bit, pos_in_guide, d, bD, bR, true, bulType, t,
					  inGuide, targetOfGuide, vTOG, vIG, bT, dirs, ind, mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);

	if (ti == (len_guide + bulDNA - bD))
	{
		saveIndices(inGuide_bit, targetOfGuide_bit, p, pos_tree, d, bD, bR, bulType, t, inGuide, targetOfGuide, vTOG, vIG, bT, dirs, ind, mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);
	}
	else if (p[pos_tree].eqkid > 0)
	{ // go to eqkid

		inGuide[gi] = *s; // save guide character
		inGuide_bit[gi] = guide_bit[pos_in_guide];
		gi++;
		if ((guide_bit[pos_in_guide] & p[pos_tree].splitchar_bit) != 0)
		{															 //MATCH CASE
			targetOfGuide[ti] = p[pos_tree].splitchar; // save target character uppercase
			targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit;
			ti++;
			nearsearch(s + 1, ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, guide_bit, pos_in_guide + 1, d, bD, bR, true, bulType, t,
						  inGuide, targetOfGuide, vTOG, vIG, bT, dirs, ind, mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);
			ti--;
		}
		else if (d > 0)
		{																	// MISMATCH case
			targetOfGuide[ti] = p[pos_tree].splitchar + 32; // save target character lowercase
			targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit;
			ti++;
			nearsearch(s + 1, ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, guide_bit, pos_in_guide + 1, d - 1, bD, bR, true, bulType, t,
						  inGuide, targetOfGuide, vTOG, vIG, bT, dirs, ind, mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);
			ti--;
		}

		if (bR > 0 && (pos_in_guide + 1) < len_guide && bulType < 1)
		{															 // BULGE RNA case
			targetOfGuide[ti] = '-';						 // update last target character with '-'
			targetOfGuide_bit[ti] = bitset<4>("0000"); // char '-'
			ti++;
			nearsearch(s + 1, ti, gi, inGuide_bit, targetOfGuide_bit, p, pos_tree, guide_bit, pos_in_guide + 1, d, bD, bR - 1, false, bulType - 1, t,
						  inGuide, targetOfGuide, vTOG, vIG, bT, dirs, ind, mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);

			ti--;
		}
		if (bD > 0 && bulType > -1)
		{															 // BULGE DNA case
			inGuide[gi - 1] = '-';							 // update last guide character with '-'
			inGuide_bit[gi - 1] = bitset<4>("0000");	//char '-'  //i've alredy inserted a char into inGuide, so it's updated with '-'
			targetOfGuide[ti] = p[pos_tree].splitchar; // save target character uppercase
			targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit;
			ti++;
			nearsearch(s, ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, guide_bit, pos_in_guide, d, bD - 1, bR, true, bulType + 1, t,
						  inGuide, targetOfGuide, vTOG, vIG, bT, dirs, ind, mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);
			ti--;
		}
		gi--;
	}
	else if (p[pos_tree].eqkid < 1 && pos_in_guide == (len_guide - 1))
	{ //current node is a leaf

		inGuide[gi] = *s; // save guide character
		inGuide_bit[gi] = guide_bit[pos_in_guide];
		gi++;
		targetOfGuide[ti] = *s == p[pos_tree].splitchar ? p[pos_tree].splitchar : p[pos_tree].splitchar + 32; // save target character

		targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit;
		ti++;
		d = (guide_bit[pos_in_guide] & p[pos_tree].splitchar_bit) != 0 ? d : d - 1; // update distance
		if (d > -1 && ti == (len_guide + bulDNA - bD))
		{ // save results
			saveIndices(inGuide_bit, targetOfGuide_bit, p, pos_tree, d, bD, bR, bulType, t, inGuide, targetOfGuide, vTOG, vIG, bT, dirs, ind, mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);
		}
		d = (guide_bit[pos_in_guide] & p[pos_tree].splitchar_bit) != 0 ? d : d + 1;
		ti--;
		gi--;
	}
}

bool compareBySize(const file_stats_struct &a, file_stats_struct &b)
{
	return a.file_dim < b.file_dim;
}

int main(int argc, char **argv)
{
	DIR *d;
	dirent *dir;
	struct stat file_stats;
	vector<file_stats_struct> file_stats_vec;
	string line;
	double start, end, globalstart, globalend; // start and end time, global start and end time

	char *genome_dir = argv[1];
	ifstream fileGuide(argv[2]);
	mm = atoi(argv[3]);
	bulDNA = atoi(argv[4]);
	bulRNA = atoi(argv[5]);
	ifstream pamfile(argv[6]);
	string name_result = argv[7];
	string type_output = argv[8];
	string num_thr_s = argv[9];
	int num_thr = stoi(num_thr_s);
	//string pam_at_start_s = argv[10];										// pam is at beginning or end of guide
	pam_at_start = false;
	// if (pam_at_start_s.compare("True") == 0)
	// 	pam_at_start = true;

	if (omp_get_max_threads() < num_thr)
		num_thr = omp_get_max_threads();

	cout << "USED THREADS " << num_thr << endl;
	getline(pamfile, line);
	transform(line.begin(), line.end(), line.begin(), ::toupper); // uppercase of the pam
	int delimiter = line.find(" ");
	string pam = line.substr(0, delimiter);
	int pamlimit = stoi(line.substr(delimiter, line.length() - 1)); //number that identifies the PAM length: NNNNNNNNNNNNNNNNNNNNNGG (3)
	if (pamlimit < 0)
	{
		pam_at_start = true;
		pamlimit = pamlimit * -1;
	}
	int pamlen = pam.length(); //length of the total PAM: (NNNNNNNNNNNNNNNNNNNNNGG) is 23
	if (!pam_at_start)
	{
		pamRNA = pam.substr(pamlen - pamlimit, pamlimit);
	}
	else
	{
		pamRNA = pam.substr(0, pamlimit); // if pam_at_start is set, then PAM = TTTNNNNNNNNNNNNNNNNNNNNN 4, i select the first 4 chars
	}
	globalstart = omp_get_wtime(); // start global time
	int numGuide;
	string iguide;
	int en;
	while (getline(fileGuide, line))
	{
		transform(line.begin(), line.end(), line.begin(), ::toupper); // toUpperCase

		if (!pam_at_start)
		{
			iguide = line.substr(0, (line.size() - 1) - pamlimit); // retrive Guide
		}
		else
		{
			iguide = line.substr(pamlimit, (line.size() - 1) - pamlimit);
		}
		guideRNA_s.push_back(iguide);
		reverse(iguide.begin(), iguide.end());
		guideRNA.push_back((char *)malloc(21 * sizeof(char)));
		copy(iguide.begin(), iguide.end(), guideRNA[numGuide]); // save Guide
		guideRNA[numGuide][pamlen - pamlimit + 1] = '\0';
		numGuide++;
	}

	cout << "guide: ";
	for (int i = 0; i < guideRNA.size(); i++)
	{
		cout << guideRNA[i];
	}
	cout << endl;
	//Transform loaded guides into bitset
	for (int i = 0; i < guideRNA.size(); i++)
	{
		vector<bitset<4>> tmp(pamlen - pamlimit);
		for (int j = 0; j < (pamlen - pamlimit); j++)
		{
			tmp[j] = bitset<4>(convertCharToBitset(guideRNA[i][j]));
		}
		guideRNA_bit.push_back(tmp);
	}

	//Get all files from a directory
	d = opendir(genome_dir);
	if (d)
	{
		while ((dir = readdir(d)) != NULL)
		{
			if (!(strcmp(".", dir->d_name)))
			{
				continue;
			}
			if (!(strcmp("..", dir->d_name)))
			{
				continue;
			}
			/*
                        if (!(strcmp(".tmp", dir->d_name)))
                        {
                                continue;
                        }
                        */
			string tmp1(genome_dir);
			string tmp2(dir->d_name);
			tmp1 += "/" + tmp2;
			stat(tmp1.c_str(), &file_stats);

			file_stats_vec.push_back(file_stats_struct());
			file_stats_vec[file_stats_vec.size() - 1].file_name = dir->d_name;
			file_stats_vec[file_stats_vec.size() - 1].file_dim = (unsigned int)file_stats.st_size;
		}
	}
	closedir(d);

	sort(file_stats_vec.begin(), file_stats_vec.end(), compareBySize);

	//file to save results and profiling
	ofstream fileResults;
	ofstream fileprofiling;
	ofstream file_ext_profiling;
	ofstream file_profiling_dna;
	ofstream file_profiling_rna;

	//Check the type_output
	create_target = true;
	create_profile = true;

	string resultwriting = "r";
	string profilewriting = "p";
	if (argv[8] == resultwriting)
	{
		fileResults.open(name_result + ".targets.txt", std::ios_base::out); //out
		fileResults << "#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge Size\n";
		create_profile = false;
	}
	else if (argv[8] == profilewriting)
	{
		fileprofiling.open(name_result + ".profile.xls", std::ios_base::out);
		file_ext_profiling.open(name_result + ".extended_profile.xls", std::ios_base::out);
		file_profiling_dna.open(name_result + ".profile_dna.xls", std::ios_base::out);
		file_profiling_rna.open(name_result + ".profile_rna.xls", std::ios_base::out);

		create_target = false;
	}
	else
	{
		fileResults.open(name_result + ".targets.txt", std::ios_base::out); //out
		fileprofiling.open(name_result + ".profile.xls", std::ios_base::out);
		file_ext_profiling.open(name_result + ".extended_profile.xls", std::ios_base::out);
		file_profiling_dna.open(name_result + ".profile_dna.xls", std::ios_base::out);
		file_profiling_rna.open(name_result + ".profile_rna.xls", std::ios_base::out);
		fileResults << "#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge Size\n";
	}

	//Create matrix for profiling //num_thr = number of layers of the matrix. Each thr update his own layer
	vector<vector<vector<int>>> profiling(pamlen - pamlimit + mm + 1, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));
	//en = guide length (BP cells), mm = number mismatches allowed (1mm, 2mm ... cells), +1 = cell for 0 mm

	//Create matrix for extended profiling
	vector<vector<vector<vector<vector<int>>>>> ext_profiling(numGuide, vector<vector<vector<vector<int>>>>(mm + 1, vector<vector<vector<int>>>(4, vector<vector<int>>(en, vector<int>(num_thr, 0)))));
	//1 dim = id of guide, 2 dim = select matrix with an x number of mms in the target, 3 dim = select nucleotide (acgt), 4 dim = select position

	//Matrix for dna profiling
	vector<vector<vector<int>>> profiling_dna_mm(pamlen - pamlimit + bulDNA + (mm + 1) * 2, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));
	vector<vector<vector<int>>> profiling_dna(pamlen - pamlimit + bulDNA + (mm + 1) * 2, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));
	//additional mm + 1 is for the column (x mismatches, 1 bulge -- x mismatches, 2 bulge)

	//Matrix for rna profiling
	vector<vector<vector<int>>> profiling_rna_mm(pamlen - pamlimit + mm + 1 + mm + 1, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));
	vector<vector<vector<int>>> profiling_rna(pamlen - pamlimit + mm + 1 + mm + 1, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));

	len_guide = pamlen - pamlimit;

	int numNodes;
	int numLeaves;
	Tleaf *targetOnDNA;

	vector<string> vecTargetOfGuide, vecInGuide, bulgeType;
	vector<char> directions;
	vector<int> indices, mismatches, bulgeSize;

	//string inGuide;
	//string targetOfGuide;
	char inGuide[50];
	char targetOfGuide[50];
	int file;
	int jk;
	int i;

	int number_of_trees = file_stats_vec.size();
	int counter = 0;
	int ten_pc = ceil(number_of_trees / (10.0));
	//#PARALLEL SECTION
	//FOR ALL THE CHR
#pragma omp parallel private(numNodes, numLeaves, targetOnDNA, vecTargetOfGuide, vecInGuide, bulgeType, directions, indices, mismatches, bulgeSize, inGuide, targetOfGuide, file, jk, i) num_threads(num_thr)
	{
#pragma omp for schedule(static, 1) nowait
		for (file = 0; file < file_stats_vec.size(); file++)
		{

			//inGuide.resize(pamlen + bulDNA);
			//targetOfGuide.resize(pamlen + bulDNA);
			string gen_dir(genome_dir);
			int last_under = file_stats_vec[file].file_name.find_last_of("_");
			int first_under = file_stats_vec[file].file_name.find_first_of("_");

			string chrName = file_stats_vec[file].file_name.substr(first_under + 1, last_under - first_under - 1);

			ifstream fileTree;
			// create tree
			vector<Tnode> albero;

			//load tst
			double load_start = omp_get_wtime();
			targetOnDNA = loadTST(gen_dir + "/" + file_stats_vec[file].file_name, albero, fileTree, numNodes, numLeaves);
			double load_end = omp_get_wtime();
			cout << "load tst " << load_end - load_start << endl;

			double guidecut = 50;
			//cout << "todna: " <<targetOnDNA[2].guideDNA << endl;
			double newguide = numGuide;
			int group_guide = ceil(newguide / guidecut);

			double guide_start = omp_get_wtime();
			for (jk = 0; jk < group_guide; jk++)
			{

				int inizio = jk * guidecut;
				int fine = MIN(((jk + 1) * guidecut), numGuide);
				int mm2 = mm;
				int bulDNA2 = bulDNA;
				int bulRNA2 = bulRNA;

				vector<bitset<4>> inGuide_bit(len_guide + pamRNA.size() + bulDNA);
				vector<bitset<4>> targetOfGuide_bit(len_guide + pamRNA.size() + bulDNA);
				for (i = inizio; i < fine; i++)
				{
					int ti = 0;
					int gi = 0;
					{
						{
							nearsearch(guideRNA[i], ti, gi, inGuide_bit, targetOfGuide_bit, albero, 0, guideRNA_bit[i], 0, mm2, bulDNA2, bulRNA2, true, 0,
										  targetOnDNA, inGuide, targetOfGuide, vecTargetOfGuide, vecInGuide, bulgeType, directions, indices, mismatches, bulgeSize,
										  i, profiling, ext_profiling, profiling_dna, profiling_dna_mm, profiling_rna, profiling_rna_mm);
						}
					}
				}
				if (create_target && indices.size() > 0)
				{
#pragma omp critical
					{

						bulgeType_glb.insert(bulgeType_glb.end(), bulgeType.begin(), bulgeType.end());
						vecInGuide_glb.insert(vecInGuide_glb.end(), vecInGuide.begin(), vecInGuide.end());
						vecTargetOfGuide_glb.insert(vecTargetOfGuide_glb.end(), vecTargetOfGuide.begin(), vecTargetOfGuide.end());

						indices_glb.insert(indices_glb.end(), indices.begin(), indices.end());
						directions_glb.insert(directions_glb.end(), directions.begin(), directions.end());
						mismatches_glb.insert(mismatches_glb.end(), mismatches.begin(), mismatches.end());
						bulgeSize_glb.insert(bulgeSize_glb.end(), bulgeSize.begin(), bulgeSize.end());

						for (i = 0; i < indices.size(); i++)
						{
							chrName_glb.push_back(chrName);
						}
					}
				}

				bulgeType.clear();
				vecInGuide.clear();
				vecTargetOfGuide.clear();
				indices.clear();
				directions.clear();
				mismatches.clear();
				bulgeSize.clear();
				albero.clear();
			}

			double guide_end = omp_get_wtime();
			cout << "nearsearch " << guide_end - guide_start << endl;

			//Free memory
			for (i = 0; i < numLeaves; i++)
			{
				delete[] targetOnDNA[i].guideDNA_bit;
				delete[] targetOnDNA[i].guideDNA;
			}

			//free targetonDNA
			free(targetOnDNA);

#pragma omp atomic
			++counter;
#pragma omp critical
			{
				cout << "ANALYZING CHROMOSOME " << chrName << " (Total progress: " << fixed << std::setprecision(1) << (100.0 * counter / number_of_trees) << "%)\n";
			}
		}
	}

	//Save profiling results
	if (create_profile)
	{
		//Cols labels for profiling output file
		fileprofiling << "GUIDE\t";
		for (int i = 0; i < len_guide; i++)
		{
			fileprofiling << "BP\t";
		}
		fileprofiling << "\t"
						  << "ONT\tOFFT\tMM/OFFT\t\t";
		for (int i = 0; i <= mm; i++)
		{
			fileprofiling << i << "MM\t";
		}
		fileprofiling << "\n";

		//Cols labels for DNA profiling output file
		file_profiling_dna << "GUIDE\t";
		for (int i = 0; i < (len_guide + bulDNA); i++)
		{
			file_profiling_dna << "BP\t";
		}
		file_profiling_dna << "\t"
								 << "ONT\tOFFT\tMM/OFFT\t\t";
		for (int i = 0; i <= mm; i++)
		{
			file_profiling_dna << i << "MM(1)\t";
			file_profiling_dna << i << "MM(2)\t";
		}
		file_profiling_dna << "\n";

		//Cols labels for RNA profiling output file
		file_profiling_rna << "GUIDE\t";
		for (int i = 0; i < (len_guide); i++)
		{
			file_profiling_rna << "BP\t";
		}
		file_profiling_rna << "\t"
								 << "ONT\tOFFT\tMM/OFFT\t\t";
		for (int i = 0; i <= mm; i++)
		{
			file_profiling_rna << i << "MM(1)\t";
			file_profiling_rna << i << "MM(2)\t";
		}
		file_profiling_rna << "\n";
		for (int i = 0; i < numGuide; i++)
		{
			saveProfileGuide(guideRNA_s[i], i, mm, len_guide, bulDNA, profiling, ext_profiling, profiling_dna, profiling_dna_mm, profiling_rna, profiling_rna_mm,
								  fileprofiling, file_ext_profiling, file_profiling_dna, file_profiling_rna, num_thr);
		}
		fileprofiling.close();
		file_ext_profiling.close();
		file_profiling_dna.close();
		file_profiling_rna.close();
	}

	if (create_target)
	{
		for (int i = 0; i < indices_glb.size(); i++)
		{
			fileResults << "\t" << bulgeType_glb[i] << "\t" << vecInGuide_glb[i] << "\t" << vecTargetOfGuide_glb[i] << "\t" << chrName_glb[i] << "\t" << indices_glb[i] << "\t" << directions_glb[i] << "\t" << mismatches_glb[i] << "\t" << bulgeSize_glb[i] << "\n";
		}
	}
	fileResults.close();
	globalend = omp_get_wtime(); // end global time
	cout << "-----------------------"
		  << "\n";
	cout << "TOTAL TIME:\t" << globalend - globalstart << "\n";

	//cout << "\nC++ end" << endl;
	return 0;
}
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
#include <iomanip>		// std::setprecision
#include <sys/stat.h> //for file dimension

#include <chrono>
using namespace std;


/***************************************
 * -18/02/2020
 * 	-Patched PAM at beginning
 * 	-Saving results every 50 guides: ACTIVE
 * 	-Added column Cluster Position
 * 
 * Note: all databases created before 25/10/2019 are missing the length of the guide in the .bin file. This search is not compatible.
 * To avoid this issue, comment the line
 * fileTree.read((char *)&offset_guide_len, sizeof(int));
 * and instead write
 * offset_guide_len = 20;
 * (or the right size of the guide the tree was build)
 */


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

string pamRNA;
vector<vector<Tnode>> albero_glb;
vector<int> guideI;
vector<int> ti, gi;
//vector<Tleaf *> targetOnDNA;
vector<vector<Tleaf>> targetOnDNA;
vector<vector<bitset<4>>> inGuide_bit;
vector<vector<bitset<4>>> targetOfGuide_bit;
vector<char *> inGuide;		
vector<char *> targetOfGuide;

vector<char *> guideRNA;								// input guide RNA
vector<string> guideRNA_s;							//input guide RNA as strings;
vector<vector<bitset<4>>> guideRNA_bit; //input guide RNA in bitset
int len_guide;
int max_bulges = 0;
bool pam_at_start;

//Create matrix for profiling //num_thr = number of layers of the matrix. Each thr update his own layer
vector<vector<vector<int>>> profiling;

//Create matrix for extended profiling
vector<vector<vector<vector<vector<int>>>>> ext_profiling;

//Create matrix for extended profiling dna and rna
vector<vector<vector<vector<int>>>> ext_profiling_dna;
vector<vector<vector<vector<int>>>> ext_profiling_rna;

//Matrix for dna profiling
vector<vector<vector<int>>> profiling_dna_mm;
vector<vector<vector<int>>> profiling_dna;

//Matrix for rna profiling
vector<vector<vector<int>>> profiling_rna_mm;
vector<vector<vector<int>>> profiling_rna;


vector<vector<string>> vecTargetOfGuide, vecInGuide, bulgeType;
vector<vector<char>> directions;
vector<vector<int>> indices, mismatches, bulgeSize, cluster_position;
vector<vector<string>> chrName_glb;



/**
 * Function that converts the read character (8 bits) into a IUPAC nucleotide character and a bitset<4> variable. 
 * @param in 8 bits character read from the .bin file
 * @param pairNuc Variable that saves the IUPAC characters converted from the read 8 bits character
 * @param pairNuc_bit Variable that saves the characters using a 4-bits representation
 */

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
		unsigned char first = in & 0xF0; 	// keep only 4 left-most bits
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


/**
 * Function that populates recursively the nodes of the TST.
 */
void deSerialize(vector<Tnode> &albero, ifstream &fileTree, vector<bitset<4>> &pairNuc_bit, char (&pairNuc)[2], char &in, uint8_t &flag, int &pNode)
{

	int i = pNode;
	albero_glb[omp_get_thread_num()][i].splitchar_bit = pairNuc_bit[flag];
	albero_glb[omp_get_thread_num()][i].splitchar = pairNuc[flag];
	if (flag)
	{
		fileTree.get(in);
		readPair(pairNuc_bit, pairNuc, in);
		flag = 0;
	}
	else
		flag++;
	if (pairNuc[flag] != '0')
	{
		albero_glb[omp_get_thread_num()][i].lokid = ++pNode;
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
	if (pairNuc[flag] != '0')
	{
		albero_glb[omp_get_thread_num()][i].hikid = ++pNode;
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
	if (pairNuc[flag] == '_')
	{
		flag++;
		fileTree.read((char *)&(albero_glb[omp_get_thread_num()][i].eqkid), sizeof(int));
	}
	else
	{
		albero_glb[omp_get_thread_num()][i].eqkid = ++pNode;
		deSerialize(albero, fileTree, pairNuc_bit, pairNuc, in, flag, pNode);
	} // go to eqkid
}

int offset_guide_len; //Number to sum to the position column when i create a tree for a 20 len guide and i search a guide with != 20 len

/**
 * Function that recreates the TST starting from the .bin file. In the first part, the array of leaves is read and the PAM are converted
 * from a 4 bits representation to a IUPAC nucleotide. After that, the deSerialize function is called to recreate the TST
 * @param path The path of the .bin file
 * @param albero Vector that contains the nodes of the TST
 * @param fileTree Variable for reading the characters of the .bin file
 * @param numNodes Number of nodes of the TST
 * @param numLeaves Number of leaves of the TST
 */ 
void loadTST(string path, vector<Tnode> &albero, ifstream &fileTree, int &numNodes, int &numLeaves) //era Tleaf *
{
	int thr = omp_get_thread_num();
	fileTree.open(path, ios::in | ios::binary);

	fileTree.read((char *)&numLeaves, sizeof(int));					// read number of leaves
	targetOnDNA[thr].resize(numLeaves);
	
	fileTree.read((char *)&offset_guide_len, sizeof(int));
	//offset_guide_len = 20;
	offset_guide_len = offset_guide_len - len_guide;
	
	
	char in;
	for (int i = 0; i < numLeaves; i++)
	{																// fill array of targets on DNA
		fileTree.read((char *)&targetOnDNA[thr][i].guideIndex, sizeof(int));  // read index of target on DNA
	
		targetOnDNA[thr][i].guideDNA_bit  = new bitset<4>[pamRNA.size()];
		targetOnDNA[thr][i].guideDNA = new char[pamRNA.size() + 1];    		// initialize PAM size
		targetOnDNA[thr][i].guideDNA[pamRNA.size()] = '\0';
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
				targetOnDNA[thr][i].guideDNA[j] = 'A';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("0001");
				//cout << "A" ;
				break;
			case 0x20:
				targetOnDNA[thr][i].guideDNA[j] = 'C';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("0010");
				//cout << "C";
				break;
			case 0x40:
				targetOnDNA[thr][i].guideDNA[j] = 'G';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("0100");
				//cout << "G" ;
				break;
			case 0x80:
				targetOnDNA[thr][i].guideDNA[j] = 'T';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("1000");
				// cout << "T" ;
				break;
			case 0x50: //R
				targetOnDNA[thr][i].guideDNA[j] = 'R';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("0101");
				// cout << "R" ;
				break;
			case 0xA0: //Y
				targetOnDNA[thr][i].guideDNA[j] = 'Y';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("1010");
				//cout << "Y";
				break;
			case 0x60: //S
				targetOnDNA[thr][i].guideDNA[j] = 'S';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("0110");
				//cout << "S" ;
				break;
			case 0x90: //W
				targetOnDNA[thr][i].guideDNA[j] = 'W';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("1001");
				//cout << "W" ;
				break;
			case 0xC0: //K
				targetOnDNA[thr][i].guideDNA[j] = 'K';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("1100");
				//cout << "K" ;
				break;
			case 0x30: //M
				targetOnDNA[thr][i].guideDNA[j] = 'M';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("0011");
				//cout << "M" ;
				break;
			case 0xE0: //B
				targetOnDNA[thr][i].guideDNA[j] = 'B';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("1110");
				//cout << "B" ;
				break;
			case 0xD0: //D
				targetOnDNA[thr][i].guideDNA[j] = 'D';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("1101");
				//cout << "D" ;
				break;
			case 0xB0: //H
				targetOnDNA[thr][i].guideDNA[j] = 'H';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("1011");
				//cout << "H";
				break;
			case 0x70: //V
				targetOnDNA[thr][i].guideDNA[j] = 'V';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("0111");
				//cout << "V";
				break;
			default:
				cerr << "Error on reading the PAM, input char is not an IUPAC nucleotide: " << mask << endl;
				break;
			}
			k++;
		}

		fileTree.get(in); 		// read index of next PAM with same guide

		if (in == '0')
		{
			targetOnDNA[thr][i].next = 0;
		}
		else
		{
			fileTree.read((char *)&targetOnDNA[thr][i].next, sizeof(int));
		}
	}

	fileTree.read((char *)&numNodes, sizeof(int)); // read number of nodes

	vector<bitset<4>> pairNuc_bit(2);
	char pairNuc[2];
	uint8_t flag = 1;

	albero_glb[omp_get_thread_num()].resize(numNodes);
	int pNode = 0;
	fileTree.get(in);

	readPair(pairNuc_bit, pairNuc, in);

	flag = 0;
	deSerialize(albero, fileTree, pairNuc_bit, pairNuc, in, flag, pNode); // deserialize TST

	fileTree.close();
}

int bulDNA, bulRNA, mm;

bool create_target;
bool create_profile;

double save_function, save_function_start, save_function_end;




/**
 * Function that saves the found target (in the DNA). It recursively go into the leaf of the current TST branch where the target is found, and
 * creates the strings for the Guide and the Target. The appropriate array are filled with data about the Target (mismatches, bulges, position in
 * the DNA etc.), for all the occurencies of the Target in the DNA. 
 * The detailedOutput functions then provide the profiling of the Guide and Target found.
 * @param p Current node of the TST
 * @param d Number of available mismatches
 * @param bD Number of available DNA bulges
 * @param bR Number of available RNA bulges
 * @param bulType Integer value that represents the type of the bulge: <0 is a RNA bulge, >0 is a DNA bulge, =0 means there are no bulges
 * (" X ") 
 * @param inGuide Array of char that contains the Guide in input and possible bulges
 * @param targetOfGuide Array of char that contains the Target found and possible bulges
 */ 
void saveIndices( Tnode *p, int d, int bD, int bR, int bulType)		
{
	if (p->lokid > 0)
		saveIndices(  &albero_glb[omp_get_thread_num()][p->lokid], d, bD, bR, bulType ); // go to lokid
	if (p->hikid > 0)
		saveIndices( &albero_glb[omp_get_thread_num()][p->hikid], d, bD, bR, bulType); // go to hikid

	if (p->eqkid < 0)
	{ // node is a leaf, save index and strand
		int thr = omp_get_thread_num();

		string g(inGuide[thr]);
		g = g.substr(0, len_guide + bulDNA - bD);

		if (pam_at_start)
		{

			for (int i = 0; i < pamRNA.size(); i++)
				g.insert(0, "N");
		}
		else
		{
			reverse(g.begin(), g.end());
			for (int i = 0; i < pamRNA.size(); i++)
				g += "N";
		}

		vector<bitset<4>> g_bit(len_guide + bulDNA);
		int j = 0;

		for (j = 0; j < (len_guide + bulDNA - bD); j++)
		{
			g_bit[j] = inGuide_bit[thr][j];
		}
		if (!pam_at_start)
			reverse(g_bit.begin(), g_bit.begin() + len_guide + bulDNA - bD);
		
		int index = (p->eqkid + 1) * -1;
		do
		{	
			
			string t(targetOfGuide[thr]);
			
			t = t.substr(0, len_guide + bulDNA - bD);
			
			
			
			

			if (pam_at_start)
			{
				string tmp(targetOnDNA[thr][index].guideDNA);
				t.insert(0, tmp); 
									
			}
			else
			{
				reverse(t.begin(), t.end());
				t += targetOnDNA[thr][index].guideDNA;
			}
			
			
			vector<bitset<4>> t_bit(len_guide +  bulDNA);
			int i = 0;
			for (i = 0; i < (len_guide + bulDNA - bD); i++)
			{
				t_bit[i] = targetOfGuide_bit[thr][i];
			}

			if (!pam_at_start){
				reverse(t_bit.begin(), t_bit.begin() + len_guide + bulDNA - bD);
			}

			//Check what type of result to create
			if (create_target)
			{

				mismatches[thr].emplace_back(mm - d); // save mismatches
				if (bulType == 0)
				{ // NO BULGE case
					bulgeType[thr].emplace_back("X");
					bulgeSize[thr].emplace_back(0);
				}
				else
				{ // BULGE case
					bulgeType[thr].emplace_back(bulType < 0 ? "RNA" : "DNA");
					bulgeSize[thr].emplace_back(bulType < 0 ? bulRNA - bR : bulDNA - bD);
				}
				if (targetOnDNA[thr][index].guideIndex < 0)
				{ // negative strand
					indices[thr].emplace_back(targetOnDNA[thr][index].guideIndex * -1);
					cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex * -1);
					if (pam_at_start){
						directions[thr].emplace_back('+');
					}
					else
					{
						directions[thr].emplace_back('-');
					}

				}
				else
				{ // strand positive
					indices[thr].emplace_back(targetOnDNA[thr][index].guideIndex + max_bulges - (bulDNA - bD) + (bulRNA - bR) + offset_guide_len );
					if (bulType == 0){		// X case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + max_bulges - (bulDNA - bD) + (bulRNA - bR) + offset_guide_len);
					}
					else if (bulType < 0){		//RNA case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + max_bulges - (bulDNA - bD) + (bulRNA - bR) + offset_guide_len - (bulRNA - bR) );
					}
					else{					//DNA case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + max_bulges - (bulDNA - bD) + (bulRNA - bR) + offset_guide_len + (bulDNA - bD) );
					}

					if (pam_at_start){
						directions[thr].emplace_back('-');
					}
					else
					{
						directions[thr].emplace_back('+');
					}
				}

				

				//Profiling and profiling ext
				if (create_profile)
				{
					if (bulType == 0)
						detailedOutputFast(guideI[thr], g_bit, t_bit, g, t, bulType, mm, len_guide, profiling, ext_profiling, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start);
					else if (bulType > 0)
						detailedOutputFastBulgeDNA(guideI[thr], g_bit, t_bit, g, t, bulType, mm, max_bulges, len_guide, bD, bulDNA, profiling_dna, profiling_dna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling);
					else
						detailedOutputFastBulgeRNA(guideI[thr], g_bit, t_bit, g, t, bulType, mm, max_bulges, len_guide, bD, bulDNA, profiling_rna, profiling_rna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_rna, ext_profiling);
					
				}
				else
				{
					vecInGuide[thr].emplace_back(g);				// save guide
					vecTargetOfGuide[thr].emplace_back(t); 			// save target
				}
			}
			else
			{ //no target, only profile
				if (bulType == 0)
					detailedOutputFast(guideI[thr], g_bit, t_bit, g, t, bulType, mm, len_guide, profiling, ext_profiling, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start);
				else if (bulType > 0)
					detailedOutputFastBulgeDNA(guideI[thr], g_bit, t_bit, g, t, bulType, mm, max_bulges, len_guide, bD, bulDNA, profiling_dna, profiling_dna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling);
				else
					detailedOutputFastBulgeRNA(guideI[thr], g_bit, t_bit, g, t, bulType, mm, max_bulges, len_guide, bD, bulDNA, profiling_rna, profiling_rna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_rna, ext_profiling);	
			}
			
			index = (targetOnDNA[thr][index].next + 1) * -1;

		} while (index > -1);
	}
	else if (p->eqkid > 0)
		saveIndices( &albero_glb[omp_get_thread_num()][p->eqkid], d, bD, bR, bulType); //go to eqkid
}

void saveIndices2( Tnode *p, int d, int bD, int bR, int bulType)	
{
	
	if (p->eqkid < 0)
	{ // node is a leaf, save index and strand
		int thr = omp_get_thread_num();

		string g(inGuide[thr]);
		g = g.substr(0, len_guide + bulDNA - bD);

		if (pam_at_start)
		{
			for (int i = 0; i < pamRNA.size(); i++)
				g.insert(0, "N");
		}
		else
		{
			reverse(g.begin(), g.end());
			for (int i = 0; i < pamRNA.size(); i++)
				g += "N";
		}

		vector<bitset<4>> g_bit(len_guide + bulDNA);
		int j = 0;

		for (j = 0; j < (len_guide + bulDNA - bD); j++)
		{
			g_bit[j] = inGuide_bit[thr][j];
		}

		if (!pam_at_start)
			reverse(g_bit.begin(), g_bit.begin() + len_guide + bulDNA - bD);
		
		int index = (p->eqkid + 1) * -1;

		do
		{	
		
			string t(targetOfGuide[thr]);
			
			t = t.substr(0, len_guide + bulDNA - bD);
			
			

			
			if (pam_at_start)
			{
				string tmp(targetOnDNA[thr][index].guideDNA);
				t.insert(0, tmp);									
			}
			else
			{
				reverse(t.begin(), t.end());
				t += targetOnDNA[thr][index].guideDNA;
			}
			
			vector<bitset<4>> t_bit(len_guide +  bulDNA);
			int i = 0;
			for (i = 0; i < (len_guide + bulDNA - bD); i++)
			{
				t_bit[i] = targetOfGuide_bit[thr][i];
			}


			if (!pam_at_start)
				reverse(t_bit.begin(), t_bit.begin() + len_guide + bulDNA - bD);
			
			//Check what type of result to create
			if (create_target)
			{

				mismatches[thr].emplace_back(mm - d); // save mismatches
				if (bulType == 0)
				{ // NO BULGE case
					bulgeType[thr].emplace_back("X");
					bulgeSize[thr].emplace_back(0);
					//TODO inserire cluster position
				}
				else
				{ // BULGE case
					bulgeType[thr].emplace_back(bulType < 0 ? "RNA" : "DNA");
					bulgeSize[thr].emplace_back(bulType < 0 ? bulRNA - bR : bulDNA - bD);
					//TODO inserire cluster position
				}
				if (targetOnDNA[thr][index].guideIndex < 0)
				{ // negative strand
					indices[thr].emplace_back(targetOnDNA[thr][index].guideIndex * -1);
					cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex * -1);
					if (pam_at_start){
						directions[thr].emplace_back('+');
					}
					else
					{
						directions[thr].emplace_back('-');
					}
				}
				else
				{ // strand positive
					indices[thr].emplace_back(targetOnDNA[thr][index].guideIndex + max_bulges - (bulDNA - bD) + (bulRNA - bR) + offset_guide_len);
					if (bulType == 0){		// X case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + max_bulges - (bulDNA - bD) + (bulRNA - bR) + offset_guide_len);
					}
					else if (bulType < 0){		//RNA case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + max_bulges - (bulDNA - bD) + (bulRNA - bR) + offset_guide_len - (bulRNA - bR) );
					}
					else{					//DNA case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + max_bulges - (bulDNA - bD) + (bulRNA - bR) + offset_guide_len + (bulDNA - bD) );
					}

					if (pam_at_start){
						directions[thr].emplace_back('-');
					}
					else
					{
						directions[thr].emplace_back('+');
					}
				}

				

				//Profiling and profiling ext
				if (create_profile)
				{
					if (bulType == 0)
						detailedOutputFast(guideI[thr], g_bit, t_bit, g, t, bulType, mm, len_guide, profiling, ext_profiling, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start);
					else if (bulType > 0)
						detailedOutputFastBulgeDNA(guideI[thr], g_bit, t_bit, g, t, bulType, mm, max_bulges, len_guide, bD, bulDNA, profiling_dna, profiling_dna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling);
					else
						detailedOutputFastBulgeRNA(guideI[thr], g_bit, t_bit, g, t, bulType, mm, max_bulges, len_guide, bD, bulDNA, profiling_rna, profiling_rna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_rna, ext_profiling);
					
				}
				else
				{
					vecInGuide[thr].emplace_back(g);				// save guide
					vecTargetOfGuide[thr].emplace_back(t); 			// save target
				}
			}
			else
			{ //no target, only profile
				if (bulType == 0)
					detailedOutputFast(guideI[thr], g_bit, t_bit, g, t, bulType, mm, len_guide, profiling, ext_profiling, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start);
				else if (bulType > 0)
					detailedOutputFastBulgeDNA(guideI[thr], g_bit, t_bit, g, t, bulType, mm, max_bulges, len_guide, bD, bulDNA, profiling_dna, profiling_dna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling);
				else
					detailedOutputFastBulgeRNA(guideI[thr], g_bit, t_bit, g, t, bulType, mm, max_bulges, len_guide, bD, bulDNA, profiling_rna, profiling_rna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_rna, ext_profiling);	
			}
			
			index = (targetOnDNA[thr][index].next + 1) * -1;

		} while (index > -1);
	}
	else if (p->eqkid > 0)
		saveIndices2( &albero_glb[omp_get_thread_num()][p->eqkid], d, bD, bR, bulType); //go to eqkid
}





// CONTROLLARE LA QUESTIONE BULGE CONTEMPORANEI IN DNA ED RNA
/**
 * Function that, given a Guide in input, searches it inside the TST. At first it goes recursively to the left child, then to the right child
 * and then compares the node character with the current Guide character.
 * 
 */ 
void nearsearch(char *s, Tnode *p,  int pos_in_guide, int d, int bD, int bR, bool goToLoHi, int bulType,
								  const int thr )
{ 
	if (p->lokid > 0 && goToLoHi && (d > 0 || bD > 0 || bR > 0 || *s < p->splitchar)){ // go to lokid
		nearsearch(s,    &albero_glb[thr][p->lokid],  pos_in_guide, d, bD, bR, true, bulType, thr);
	}

	if (p->hikid > 0 && goToLoHi && (d > 0 || bD > 0 || bR > 0 || *s > p->splitchar)) {// go to hikid
		nearsearch(s,    &albero_glb[thr][p->hikid],  pos_in_guide, d, bD, bR, true, bulType, thr );
		
	}

	if (ti[thr] == (len_guide + bulDNA - bD))
	{
		
		save_function_start = omp_get_wtime();
		saveIndices( p, d, bD, bR, bulType);
		save_function_end = omp_get_wtime();
		save_function += save_function_end - save_function_start;
	}
	else if (p->eqkid > 0)
	{ // go to eqkid

		inGuide[thr][gi[thr]] = *s; // save guide character
		inGuide_bit[thr][gi[thr]] = guideRNA_bit[guideI[thr]][pos_in_guide];
		gi[thr]++;
		if ((guideRNA_bit[guideI[thr]][pos_in_guide] & p->splitchar_bit) != 0)
		{																						 //MATCH CASE
			targetOfGuide[thr][ti[thr]] = p->splitchar; // save target character uppercase
			targetOfGuide_bit[thr][ti[thr]] = p->splitchar_bit;
			ti[thr]++;
			nearsearch(s + 1, &albero_glb[thr][p->eqkid],  pos_in_guide + 1, d, bD, bR, true, bulType, thr);
			ti[thr]--;
		}
		else if (d > 0)
		{																									// MISMATCH case
			targetOfGuide[thr][ti[thr]] = p->splitchar + 32; // save target character lowercase
			targetOfGuide_bit[thr][ti[thr]] = p->splitchar_bit;
			ti[thr]++;
			nearsearch(s + 1,  &albero_glb[thr][p->eqkid],  pos_in_guide + 1, d - 1, bD, bR, true, bulType, thr);
			ti[thr]--;
		}

		if (bR > 0 && (pos_in_guide + 1) < len_guide && bulType < 1)
		{																						 // BULGE RNA case
			targetOfGuide[thr][ti[thr]] = '-';									 // update last target character with '-'
			targetOfGuide_bit[thr][ti[thr]] = bitset<4>("0000"); // char '-'
			ti[thr]++;
			nearsearch(s + 1,   p,  pos_in_guide + 1, d, bD, bR - 1, false, bulType - 1, thr);
			ti[thr]--;
		}
		if (bD > 0 && bulType > -1)
		{																						 // BULGE DNA case
			inGuide[thr][gi[thr] - 1] = '-';										 // update last guide character with '-'
			inGuide_bit[thr][gi[thr] - 1] = bitset<4>("0000");	 //char '-'  //i've alredy inserted a char into inGuide, so it's updated with '-'
			targetOfGuide[thr][ti[thr]] = p->splitchar; // save target character uppercase
			targetOfGuide_bit[thr][ti[thr]] = p->splitchar_bit;
			ti[thr]++;
			nearsearch(s,   &albero_glb[thr][p->eqkid],  pos_in_guide, d, bD - 1, bR, true, bulType + 1, thr);
			
			ti[thr]--;
		}
		gi[thr]--;
	}
	else if (p->eqkid < 1 && pos_in_guide == (len_guide - 1))
	{ //current node is a leaf
		//Questa porzione di codice serve (forse) per evitare di avere un bulge a inizio guida, eg -CTAAC...NNN
		
		
		inGuide[thr][gi[thr]] = *s; // save guide character
		inGuide_bit[thr][gi[thr]] = guideRNA_bit[guideI[thr]][pos_in_guide];
		gi[thr]++;
		targetOfGuide[thr][ti[thr]]  =  (guideRNA_bit[guideI[thr]][pos_in_guide] & p->splitchar_bit) != 0 ? p->splitchar : p->splitchar +32;

		targetOfGuide_bit[thr][ti[thr]] = p->splitchar_bit;
		ti[thr]++;


		d = (guideRNA_bit[guideI[thr]][pos_in_guide] & p->splitchar_bit) != 0 ? d : d - 1; // update distance


		if (d > -1 && ti[thr] == (len_guide + bulDNA - bD))
		{ // save results
			
			
			save_function_start = omp_get_wtime(); 
			
			//saveIndices( p, d, bD, bR, bulType);
			saveIndices2(p, d, bD, bR, bulType);
			
			save_function_end = omp_get_wtime();
			save_function += save_function_end - save_function_start;
		}
		d = (guideRNA_bit[guideI[thr]][pos_in_guide] & p->splitchar_bit) != 0 ? d : d + 1;
		ti[thr]--;
		gi[thr]--;
		
	}
}

bool compareBySize(const file_stats_struct &a, file_stats_struct &b)
{
	return a.file_dim < b.file_dim;
}

int main(int argc, char **argv)
{
	// cout << "VERSIONE DI TEST PER PAM INIZIO -> done"<< endl;
	// cout << "VERSIONE DI TEST PER CONTROLLO SALVATAGGIO OGNI 50 GUIDE -> (attiva)"<< endl;
	// cout << "VERSIONE DI TEST PER AGGIUNTA COLONNA CLUSTER POSITION -> done"  << endl;
	DIR *d;
	dirent *dir;
	struct stat file_stats;
	vector<file_stats_struct> file_stats_vec;
	string line;
	double globalstart, globalend; // start and end time, global start and end time

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
	pam_at_start = false;
	
	
	string max_bulges_string = argv[10];
	max_bulges = stoi(max_bulges_string);	
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
		pamRNA = pam.substr(0, pamlimit); // if pam_at_start is set, then PAM = TTTNNNNNNNNNNNNNNNNNNNNN -4, i select the first 4 chars
	}
	globalstart = omp_get_wtime(); // start global time
	int numGuide;
	string iguide;
	while (getline(fileGuide, line))
	{
		transform(line.begin(), line.end(), line.begin(), ::toupper); // toUpperCase

		if (!pam_at_start)
		{
			iguide = line.substr(0, pamlen - pamlimit); // retrive Guide
		}
		else
		{
			iguide = line.substr(pamlimit, pamlen - pamlimit);
		}

		guideRNA_s.push_back(iguide);
		if (!pam_at_start) {
			reverse(iguide.begin(), iguide.end());
		}
		guideRNA.push_back((char *)malloc((pamlen-pamlimit) * sizeof(char)));
		copy(iguide.begin(), iguide.end(), guideRNA[numGuide]); // save Guide
		guideRNA[numGuide][pamlen - pamlimit + 1] = '\0';
		numGuide++;
	}

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
	ofstream file_profiling_complete;
	//Check the type_output
	create_target = true;
	create_profile = true;

	string resultwriting = "r";
	string profilewriting = "p";
	if (argv[8] == resultwriting)
	{
		fileResults.open(name_result + ".targets.txt", std::ios_base::out); //out
		fileResults << "#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge Size\tTotal\n";
		create_profile = false;
	}
	else if (argv[8] == profilewriting)
	{
		fileprofiling.open(name_result + ".profile.xls", std::ios_base::out);
		file_ext_profiling.open(name_result + ".extended_profile.xls", std::ios_base::out);
		file_profiling_dna.open(name_result + ".profile_dna.xls", std::ios_base::out);
		file_profiling_rna.open(name_result + ".profile_rna.xls", std::ios_base::out);
		file_profiling_complete.open(name_result + ".profile_complete.xls", std::ios_base::out);
		create_target = false;
	}
	else
	{
		fileResults.open(name_result + ".targets.txt", std::ios_base::out); //out
		fileprofiling.open(name_result + ".profile.xls", std::ios_base::out);
		file_ext_profiling.open(name_result + ".extended_profile.xls", std::ios_base::out);
		file_profiling_dna.open(name_result + ".profile_dna.xls", std::ios_base::out);
		file_profiling_rna.open(name_result + ".profile_rna.xls", std::ios_base::out);
		file_profiling_complete.open(name_result + ".profile_complete.xls", std::ios_base::out);
		fileResults << "#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge Size\tTotal\n";
	}

	//Resize matrix for profiling //num_thr = number of layers of the matrix. Each thr update his own layer
	profiling.resize(pamlen - pamlimit + mm + 1, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));
	//en = guide length (BP cells), mm = number mismatches allowed (1mm, 2mm ... cells), +1 = cell for 0 mm

	//Create matrix for extended profiling
	ext_profiling.resize(numGuide, vector<vector<vector<vector<int>>>>(mm + 1, vector<vector<vector<int>>>(4, vector<vector<int>>(pamlen - pamlimit, vector<int>(num_thr, 0)))));
	//1 dim = id of guide, 2 dim = select matrix with an x number of mms in the target, 3 dim = select nucleotide (acgt), 4 dim = select position

	//Resize matrix for extended profiling dna and rna
	ext_profiling_dna.resize(numGuide, vector<vector<vector<int>>>(mm + 1, vector<vector<int>>(pamlen - pamlimit, vector<int>(num_thr, 0))));
	ext_profiling_rna.resize(numGuide, vector<vector<vector<int>>>(mm + 1, vector<vector<int>>(pamlen - pamlimit, vector<int>(num_thr, 0))));
	
	//Matrix for dna profiling
	profiling_dna_mm.resize(pamlen - pamlimit + bulDNA + (mm + 1) * max_bulges, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));
	profiling_dna.resize(pamlen - pamlimit + bulDNA + (mm + 1) * max_bulges, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));
	//additional mm + 1 is for the column (x mismatches, 1 bulge -- x mismatches, 2 bulge)

	//Matrix for rna profiling
	profiling_rna_mm.resize(pamlen - pamlimit + (mm + 1) * max_bulges , vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));
	profiling_rna.resize(pamlen - pamlimit + (mm + 1)* max_bulges, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));


	vecTargetOfGuide.resize(num_thr);
	vecInGuide.resize(num_thr); 
	bulgeType.resize(num_thr);
	directions.resize(num_thr);
	indices.resize(num_thr); 
	cluster_position.resize(num_thr);
	mismatches.resize(num_thr);
	bulgeSize.resize(num_thr);
	chrName_glb.resize(num_thr);
	guideI.resize(num_thr);
	targetOnDNA.resize(num_thr);
	albero_glb.resize(num_thr);
	ti.resize(num_thr);
	gi.resize(num_thr);
	inGuide_bit.resize(num_thr);
	targetOfGuide_bit.resize(num_thr);
	len_guide = pamlen - pamlimit;
	inGuide.resize(num_thr);
	targetOfGuide.resize(num_thr);
	
	int numNodes;
	int numLeaves;

	char c_inGuide[50];
	char c_targetOfGuide[50];
	int file;
	int jk;
	int i;

	int number_of_trees = file_stats_vec.size();
	int counter = 0;
	
	//#PARALLEL SECTION
	//FOR ALL THE CHR
#pragma omp parallel private(c_inGuide, c_targetOfGuide, numNodes, numLeaves, file, jk, i) num_threads(num_thr)
	{
#pragma omp for schedule(guided) nowait
		for (file = 0; file < file_stats_vec.size(); file++)
		{
			save_function=0;
			numNodes = 0;
			int thr = omp_get_thread_num();
			inGuide[thr] = c_inGuide;
			targetOfGuide[thr] = c_targetOfGuide;

			string gen_dir(genome_dir);
			int last_under = file_stats_vec[file].file_name.find_last_of("_");
			int first_under = file_stats_vec[file].file_name.find_first_of("_");

			string chrName = file_stats_vec[file].file_name.substr(first_under + 1, last_under - first_under - 1);
			ifstream fileTree;

			// create tree
			vector<Tnode> albero;
			//load tst
			loadTST(gen_dir + "/" + file_stats_vec[file].file_name, albero, fileTree, numNodes, numLeaves);

		
			double guidecut = 200;//50;
			double newguide = numGuide;
			int group_guide = ceil(newguide / guidecut);
			for (jk = 0; jk < group_guide; jk++)
			{

				int inizio = jk * guidecut;
				int fine = MIN(((jk + 1) * guidecut), numGuide);
				int mm2 = mm;
				int bulDNA2 = bulDNA;
				int bulRNA2 = bulRNA;

				inGuide_bit[thr].resize(len_guide + pamRNA.size() + bulDNA);
				targetOfGuide_bit[thr].resize(len_guide + pamRNA.size() + bulDNA);
				for (i = inizio; i < fine; i++)
				{
					ti[thr] = 0;
					gi[thr] = 0;
					guideI[thr] = i;
					
					nearsearch(guideRNA[i],  &albero_glb[thr][0], 0, mm2, bulDNA2, bulRNA2, true, 0, thr);					
				}

				if (create_target)
				{
					for (i = 0; i < indices[omp_get_thread_num()].size(); i++)
					{
						chrName_glb[omp_get_thread_num()].push_back(chrName);
					}
					
				}

			//}   //-> comment/uncomment to deactivate/activate saving every 50 guides
			
			//Save results
			
			if (create_target)
			{
	
				#pragma omp critical
				{
					for (int j = 0; j < indices[thr].size(); j++)
					{
						
						fileResults << bulgeType[thr][j] << "\t" << vecInGuide[thr][j] << "\t" << vecTargetOfGuide[thr][j] << "\t" << chrName_glb[thr][j] << "\t" << indices[thr][j] << "\t" << cluster_position[thr][j] << "\t" << directions[thr][j] << "\t" << mismatches[thr][j] << "\t" << bulgeSize[thr][j] << "\t" <<  mismatches[thr][j] + bulgeSize[thr][j] << "\n";

					}
				}

				bulgeType[thr].clear(); //bulgeType.resize(num_thr);
				vecInGuide[thr].clear(); //vecInGuide.resize(num_thr);
				vecTargetOfGuide[thr].clear(); //vecTargetOfGuide.resize(num_thr);
				chrName_glb[thr].clear(); //chrName_glb.resize(num_thr);
				indices[thr].clear(); //indices.resize(num_thr);
				cluster_position[thr].clear();
				directions[thr].clear(); //directions.resize(num_thr);
				mismatches[thr].clear(); //mismatches.resize(num_thr);
				bulgeSize[thr].clear(); //bulgeSize.resize(num_thr);
			}
			}			//uncomment/comment to activate/deactivate saving every 50 guides
			
			//Free memory targetOnDNA
			for (i = 0; i < numLeaves; i++)
			{
				delete[] targetOnDNA[thr][i].guideDNA_bit;
				delete[] targetOnDNA[thr][i].guideDNA;
			}

			albero_glb[thr].clear(); 
			targetOnDNA[thr].clear();

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
		for (int i = 0; i < (len_guide); i++) //+bulDNA
		{
			file_profiling_dna << "BP\t";
		}
		file_profiling_dna << "\t"
											 << "ONT\tOFFT\tMM/OFFT\t\t";
		for (int i = 0; i <= mm; i++)
		{
			for (int j = 1; j <=max_bulges; j++){
				file_profiling_dna << i << "MM(" << j << ")\t";
			}
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
			for (int j = 1; j <=max_bulges; j++){
				file_profiling_rna << i << "MM(" << j << ")\t";
			}
		}
		file_profiling_rna << "\n";

		//Cols label for file_profiling_complete
		file_profiling_complete << "GUIDE\t";
		for (int i = 0; i < (len_guide); i++)
		{
			file_profiling_complete << "BP\t";
		}
		file_profiling_complete << "\t"
											 << "ONT\tOFFT\tMM/OFFT\t\t";
		for (int i = 0; i <= mm; i++)
		{
			file_profiling_complete << i << "MM\t";
		}
		file_profiling_complete << "\n";


		for (int i = 0; i < numGuide; i++)
		{
			saveProfileGuide(guideRNA_s[i], i, mm, max_bulges, len_guide, bulDNA, profiling, ext_profiling, profiling_dna, profiling_dna_mm, profiling_rna, profiling_rna_mm,
											 fileprofiling, file_ext_profiling, file_profiling_dna, file_profiling_rna, file_profiling_complete , num_thr, pamRNA.size(), pam_at_start, ext_profiling_dna, ext_profiling_rna);
		}
		fileprofiling.close();
		file_ext_profiling.close();
		file_profiling_dna.close();
		file_profiling_rna.close();
		file_profiling_complete.close();
	}

	fileResults.close();
	globalend = omp_get_wtime(); // end global time
	cout << "-----------------------"
			 << "\n";
	cout << "TOTAL TIME:\t" << globalend - globalstart << "\n";

	return 0;
}
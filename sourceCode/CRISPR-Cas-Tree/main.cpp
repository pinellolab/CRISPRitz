#include "pamGenerator.cpp"
#include <iostream>
#include <ostream>
#include <omp.h>
#include <string>
#include <vector>
#include <fstream>
//#include <parallel/algorithm>
#include <algorithm>
#include <algorithm>
#include <unistd.h>

/***************************************
 * -18/02/2020
 * 	-Patched PAM at beginning
 * 	-Current Version: NO __gnu_parallel::
 * Note: all databases created before 25/10/2019 are missing the length of the guide in the .bin file
 * 
 */

using namespace std;

typedef struct tnode *Tptr;
typedef struct tnode
{
	char splitchar;
	int lokid, eqkid, hikid;
} Tnode;

typedef struct tleaf
{
	int guideIndex;
	string guideDNA;
	const char * guideDNA_char;
	string pamDNA;
	int next;
} Tleaf;

#define MAXWORDS 32000000
#define MAXCHARS (MAXWORDS * 30)

const int TARG_IN_GROUP = 5000000; //number of strings for each TST; Change to have smaller or bigger TSTs
					   
Tptr tree; 
int nodeUsed; // number of nodes

vector<Tleaf> targetOnDNA;		// array of target on DNA
string pamRNA;				 	// input pam RNA
string chrName;				 	// input chromosome name
string chrSeq;				 	// input chromosome sequence
int gruppo;
bool pam_at_start;


// function that inserts a guide DNA in TST
void insert(string st, int i, int i2)
{
	int d;					 		// distance between characters
	const char *s = st.c_str();
	Tptr pp = tree; 				// insert start from root
	while (nodeUsed)
	{ // move through TST
		if ((d = *s - pp->splitchar) == 0)
		{
			if (*++s == 0)
			{ // update leaf
				targetOnDNA[i].next = pp->eqkid;
				pp->eqkid = (i2 + 1) * -1;
				return;
			}
			if (!pp->eqkid)
			{ // go to eqkid
				pp->eqkid = nodeUsed;
				break;
			}
			pp = &tree[pp->eqkid];
		}
		else if (d < 0)
		{ // go to lowkid
			if (!pp->lokid)
			{
				pp->lokid = nodeUsed;
				break;
			}
			pp = &tree[pp->lokid];
		}
		else
		{
			if (!pp->hikid)
			{ // go to hikid
				pp->hikid = nodeUsed;
				break;
			}
			pp = &tree[pp->hikid];
		}
	}
	for (;;)
	{ // insert new node in TST
		pp = &tree[nodeUsed];
		pp->splitchar = *s;
		pp->lokid = pp->eqkid = pp->hikid = 0;
		nodeUsed++;
		if (*++s == 0)
		{ // add leaf information
			pp->eqkid = (i2 + 1) * -1;
			return;
		}
		pp->eqkid = nodeUsed;
	}
}

ofstream fileTree;

char pairNuc[2];
bool flag = true;
unsigned char bitNuc;
double num_iupac = 0;

//check http://cactus.io/resources/toolbox/decimal-binary-octal-hexadecimal-conversion
void writePair()
{

	if (pairNuc[0] == '_')
		bitNuc = 0xF0; //1 1 1 1   0 0 0 0
	else
	{

		switch (pairNuc[0])
		{ //fill 4 most significant bits (left-most)
		case 'A':
			bitNuc = 0x10;
			break; //0 0 0 1   0 0 0 0
		case 'C':
			bitNuc = 0x20;
			break; //0 0 1 0   0 0 0 0
		case 'G':
			bitNuc = 0x40;
			break; //0 1 0 0   0 0 0 0
		case 'T':
			bitNuc = 0x80;
			break; //1 0 0 0   0 0 0 0
		//case 'N': 						//not included because no N can appear in the targets
		case 'R':
			bitNuc = 0x50;
			num_iupac++;
			break; //0 1 0 1   0 0 0 0
		case 'Y':
			bitNuc = 0xA0;
			num_iupac++;
			break; //1 0 1 0   0 0 0 0
		case 'S':
			bitNuc = 0x60;
			num_iupac++;
			break; //0 1 1 0   0 0 0 0
		case 'W':
			bitNuc = 0x90;
			num_iupac++;
			break; //1 0 0 1   0 0 0 0
		case 'K':
			bitNuc = 0xC0;
			num_iupac++;
			break; //1 1 0 0   0 0 0 0
		case 'M':
			bitNuc = 0x30;
			num_iupac++;
			break; //0 0 1 1   0 0 0 0
		case 'B':
			bitNuc = 0xE0;
			num_iupac++;
			break; //1 1 1 0   0 0 0 0
		case 'D':
			bitNuc = 0xD0;
			num_iupac++;
			break; //1 1 0 1   0 0 0 0
		case 'H':
			bitNuc = 0xB0;
			num_iupac++;
			break; //1 0 1 1   0 0 0 0
		case 'V':
			bitNuc = 0x70;
			num_iupac++;
			break; //0 1 1 1   0 0 0 0
		case '0':
			bitNuc = 0x0;
			break; //0 0 0 0   0 0 0 0
		default:
			cerr << "The character (" << pairNuc[0] << ") is not part of the IUPAC nucleotide nomenclature" << endl;
			bitNuc = 0x0;
			break;
		}
		switch (pairNuc[1])
		{ //fill 4 least significant bits (right-most)
		case 'A':
			bitNuc += 0x1;
			break; //0 0 0 0   0 0 0 1
		case 'C':
			bitNuc += 0x2;
			break; //0 0 0 0   0 0 1 0
		case 'G':
			bitNuc += 0x4;
			break; //0 0 0 0   0 1 0 0
		case 'T':
			bitNuc += 0x8;
			break; //0 0 0 0   1 0 0 0
		//case 'N': 						//not included because no N can appear in the targets
		case '_':
			bitNuc += 0x0F;
			break; //0 0 0 0   1 1 1 1
		case 'R':
			bitNuc += 0x05;
			num_iupac++;
			break; //0 0 0 0   0 1 0 1
		case 'Y':
			bitNuc += 0x0A;
			num_iupac++;
			break; //0 0 0 0   1 0 1 0
		case 'S':
			bitNuc += 0x06;
			num_iupac++;
			break; //0 0 0 0   0 1 1 0
		case 'W':
			bitNuc += 0x09;
			num_iupac++;
			break; //0 0 0 0   1 0 0 1
		case 'K':
			bitNuc += 0x0C;
			num_iupac++;
			break; //0 0 0 0   1 1 0 0
		case 'M':
			bitNuc += 0x03;
			num_iupac++;
			break; //0 0 0 0   0 0 1 1
		case 'B':
			bitNuc += 0x0E;
			num_iupac++;
			break; //0 0 0 0   1 1 1 0
		case 'D':
			bitNuc += 0x0D;
			num_iupac++;
			break; //0 0 0 0   1 1 0 1
		case 'H':
			bitNuc += 0x0B;
			num_iupac++;
			break; //0 0 0 0   1 0 1 1
		case 'V':
			bitNuc += 0x07;
			num_iupac++;
			break; //0 0 0 0   0 1 1 1
		case '0':
			bitNuc += 0x0;
			break; //0 0 0 0   0 0 0 0
		default:
			cerr << "The character (" << pairNuc[1] << ") is not part of the IUPAC nucleotide nomenclature" << endl;
			bitNuc += 0x0;
			break;
		}
	}
	fileTree.put(bitNuc);
}

// Serialize TST
void serialize(Tptr p)
{

	if (flag)
	{ // write current node and recur for its children
		pairNuc[0] = p->splitchar;
		flag = false;
	}
	else
	{
		pairNuc[1] = p->splitchar;
		writePair();
		flag = true;
	}
	if (p->lokid > 0)
		serialize(&tree[p->lokid]); // go to lokid
	else
	{
		if (flag)
		{
			pairNuc[0] = '0';
			flag = false;
		}
		else
		{
			pairNuc[1] = '0';
			writePair();
			flag = true;
		}
	}
	if (p->hikid > 0)
		serialize(&tree[p->hikid]); // go to hikid
	else
	{
		if (flag)
		{
			pairNuc[0] = '0';
			flag = false;
		}
		else
		{
			pairNuc[1] = '0';
			writePair();
			flag = true;
		}
	}
	if (p->eqkid > 0)
		serialize(&tree[p->eqkid]); // go to eqkid
	else
	{
		if (flag)
			pairNuc[0] = '_';
		else
			pairNuc[1] = '_';
		writePair();
		flag = true;
		fileTree.write((char *)&p->eqkid, sizeof(int)); // store array PAM index
	}
}


int len_guide_used = 20;
// Write to file
void saveTST(int inizio, int fine, int part)
{
	double start, end;
	int arrayDim = fine - inizio;
	fileTree.open(pamRNA + "_" + chrName + "_" + to_string(part) + ".bin", ios::out | ios::binary);

	fileTree.write((char *)&arrayDim, sizeof(int)); // write number of targets
	fileTree.write((char*)&len_guide_used, sizeof(int)); //write len of guide etc etc
	for (int i = inizio; i < fine; i++)
	{																	 // write array of targets on DNA
		fileTree.write((char *)&targetOnDNA[i].guideIndex, sizeof(int)); // write index of target on DNA
		int k = 0;
		bitNuc = 0;
		int counter = 0;
		const char * ppp = targetOnDNA[i].pamDNA.c_str(); 
		do
		{ // write target site PAM
			counter++;
			switch (*ppp)
			{ //bits table: check the writePair() function
			case 'A':
				bitNuc += 0x1;
				break;
			case 'C':
				bitNuc += 0x2;
				break;
			case 'G':
				bitNuc += 0x4;
				break;
			case 'T':
				bitNuc += 0x8;
				break;
			case 'R':
				bitNuc += 0x05;
				break;
			case 'Y':
				bitNuc += 0x0A;
				break;
			case 'S':
				bitNuc += 0x06;
				break;
			case 'W':
				bitNuc += 0x09;
				break;
			case 'K':
				bitNuc += 0x0C;
				break;
			case 'M':
				bitNuc += 0x03;
				break;
			case 'B':
				bitNuc += 0x0E;
				break;
			case 'D':
				bitNuc += 0x0D;
				break;
			case 'H':
				bitNuc += 0x0B;
				break;
			case 'V':
				bitNuc += 0x07;
				break;
			default:
				cerr << "The PAM cannot contain the character (" << *ppp << ")" << endl;
				break;
			}
			k++;
			*ppp++;
			if (!*ppp || k == 2)
			{
				if (counter % 3 == 0 && counter != 0)
					bitNuc <<= 4;
				fileTree.put(bitNuc);
				bitNuc = 0;
				k = 0;
			}

			bitNuc <<= 4;
		} while (*ppp);

		if (targetOnDNA[i].next)
		{ // write index of next
			fileTree.put('_');
			fileTree.write((char *)&targetOnDNA[i].next, sizeof(int));
		}
		else
		{
			fileTree.put('0');
		}
	}
	fileTree.write((char *)&nodeUsed, sizeof(int)); // write number of nodes
	serialize(&tree[0]);							// serialize TST

	fileTree.close();
}

/* TIMING */

//Sorting function
bool compareFunc(Tleaf a, Tleaf b)
{
	return a.guideDNA.compare(b.guideDNA) < 0;
}

//Insert string in the tree
void insall(int l, int r)
{
	if (r < l)
		return;
	int m = l + (r - l) / 2;
	insert(targetOnDNA[m].guideDNA, m,(m - (gruppo -1) * TARG_IN_GROUP ));
	insall(l, m - 1);
	insall(m + 1, r);
}

int main(int argc, char **argv)
{
	//cout << "Versione di test per pam all'inizio -> DONE" << endl;
	string line;											  // line of fasta file
	double start, end, globalstart, globalend, globalpartial; // start and end time, global start and end time
	vector<string> all_pam;									  // vector of all possible pam RNA
	vector<int> pamIndices;									  // vector of target indices of pam RNA on DNA
	ifstream fasta(argv[1]);								  // input fasta file
	ifstream pamfile(argv[2]);								  // input pam.txt
	int max_bulges = stoi(argv[4]);							//max allowed bulges
	pam_at_start = false;
	globalstart = omp_get_wtime(); // start global time

	// ----------------------- READ INPUT FASTA ---------------------------
	//Read chromosome name
	getline(fasta, chrName);
	chrName = chrName.substr(1, chrName.length() - 1);
	cout << "Load chr:\t";
	start = omp_get_wtime();
	while (getline(fasta, line).good())
	{ // read chromosome sequence
		chrSeq += line;
	}
	//__gnu_parallel::transform(chrSeq.begin(), chrSeq.end(), chrSeq.begin(), ::toupper); // parallelized to uppercase
	transform(chrSeq.begin(), chrSeq.end(), chrSeq.begin(), ::toupper);
	end = omp_get_wtime();
	cout << end - start << "\n";

	// ------------------- GENERATE ALL POSSIBLE PAM -------------------
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
	 
	int pamlen = pam.length();										//length of the total PAM: (NNNNNNNNNNNNNNNNNNNNNGG) is 23
	
	len_guide_used = pamlen - pamlimit;
	if (!pam_at_start)
	{
		pamRNA = pam.substr(pamlen - pamlimit, pamlimit);
	}
	else
	{
		pamRNA = pam.substr(0, pamlimit); // if pam_at_start is set, then PAM = TTTNNNNNNNNNNNNNNNNNNNNN -4, i select the first 4 chars
	}
	all_pam = generatePam(pamRNA); // generate a vector of each possible input pam

	// make list of all possible pam RNA
	string list[all_pam.size()];

	for (int j = 0; j < all_pam.size(); j++)
	{
		list[j] = all_pam[j];
	}

	// ------------------- SEARCH PAM IN THE CHROMOSOME -------------------
	cout << "Search PAM:\t";
	start = omp_get_wtime();
	searchWords(pamIndices, list, all_pam.size(), chrSeq, pamlen, pamlimit, pam_at_start, max_bulges);
	all_pam.clear();
	end = omp_get_wtime();
	cout << end - start << "\n";

	// ------------------------ CREATE THE TST ----------------------------

	int discarded = 0;
	cout << "Retrieve seq:\t";
	start = omp_get_wtime();
	int i = 0;
	int counter_index = 0;
	targetOnDNA.resize(pamIndices.size());
	if (pam_at_start){
		for (i = 0; i < pamIndices.size(); i++)
		{
	
			string target;
			if (pamIndices[i] < 0) 		//String found in positive strand (PAM AT BEGINNING case)
			{
				target = chrSeq.substr(pamIndices[i] * -1, pamlen + max_bulges); //extract target + pam + 2 char for bulges from the chromosome
				
	
				if (target.find('N') != std::string::npos)		   //if 'N' is in the target, remove the target
				{
					counter_index--;
					discarded++;
				}
				else
				{
					
					string tmp_pam_str;
					tmp_pam_str = target.substr(0,pamRNA.length());
					reverse(tmp_pam_str.begin(), tmp_pam_str.end());
					
					targetOnDNA[counter_index] = (Tleaf){pamIndices[i], target.substr(pamRNA.length(), pamlen - pamRNA.length()+ max_bulges), target.substr(pamRNA.length(), pamlen - pamRNA.length()+ max_bulges).c_str(),
					 				tmp_pam_str,
									0}; //salvo l'indice del target
				
				}
			}
			else
			{
				target = chrSeq.substr((pamIndices[i]), pamlen + max_bulges);
			
				if (target.find('N') != std::string::npos)
				{
					counter_index--;
					discarded++;
				}
				else
				{
					
					string tmp;
					for (char &c : target) //complemento dei nucleotidi per pam negativa
						switch (c)
						{
						case 'A':
							tmp += 'T';
							break;
						case 'T':
							tmp += 'A';
							break;
						case 'C':
							tmp += 'G';
							break;
						case 'G':
							tmp += 'C';
							break;
						case 'R':
							tmp += 'Y';
							break;
						case 'Y':
							tmp += 'R';
							break;
						case 'S':
							tmp += 'S';
							break;
						case 'W':
							tmp += 'W';
							break;
						case 'M':
							tmp += 'K';
							break;
						case 'K':
							tmp += 'M';
							break;
						case 'H':
							tmp += 'D';
							break;
						case 'D':
							tmp += 'H';
							break;
						case 'B':
							tmp += 'V';
							break;
						case 'V':
							tmp += 'B';
							break;
						default:
							cerr << "The character (" << c << ") of the PAM is not part of the IUPAC nucleotide nomenclature" << endl;
							tmp += 'N';
							break;
						}

					reverse(tmp.begin(), tmp.end());
					
					string tmp_pam_str;
					tmp_pam_str = tmp.substr(0,pamRNA.length());
					reverse(tmp_pam_str.begin(), tmp_pam_str.end());

					targetOnDNA[counter_index] = (Tleaf){pamIndices[i], tmp.substr(pamRNA.length(), pamlen - pamRNA.length()+ max_bulges), tmp.substr(pamRNA.length(), pamlen - pamRNA.length()+ max_bulges).c_str(),
					 				tmp_pam_str, 
									 0}; //salvo l'indice del target
				}
			}
			
			counter_index++;
		}
	}
	else{			//PAM AT END
		for (i = 0; i < pamIndices.size(); i++)
		{
			string target;
			if (pamIndices[i] > 0) 
			{
				target = chrSeq.substr(pamIndices[i], pamlen + max_bulges); //extract target + pam + 2 char for bulges from the chromosome
								
				if (target.find('N') != std::string::npos)		   //if 'N' is in the target, remove the target
				{
					counter_index--;
					discarded++;
				}
				else
				{
					
					reverse(target.begin(), target.end()); //reverse per aggiungere nell'albero
					
					targetOnDNA[counter_index] = (Tleaf){pamIndices[i], target.substr(pamRNA.length()), target.substr(pamRNA.length()).c_str(),
									target.substr(0, pamRNA.length()), 0}; //salvo l'indice del target			
					
				}
			}
			else
			{
				target = chrSeq.substr((pamIndices[i]) * -1, pamlen + max_bulges);
				
				if (target.find('N') != std::string::npos)
				{
					counter_index--;
					discarded++;
				}
				else
				{
					
					string tmp;
					for (char &c : target) //complemento dei nucleotidi per pam negativa
						switch (c)
						{
						case 'A':
							tmp += 'T';
							break;
						case 'T':
							tmp += 'A';
							break;
						case 'C':
							tmp += 'G';
							break;
						case 'G':
							tmp += 'C';
							break;
						case 'R':
							tmp += 'Y';
							break;
						case 'Y':
							tmp += 'R';
							break;
						case 'S':
							tmp += 'S';
							break;
						case 'W':
							tmp += 'W';
							break;
						case 'M':
							tmp += 'K';
							break;
						case 'K':
							tmp += 'M';
							break;
						case 'H':
							tmp += 'D';
							break;
						case 'D':
							tmp += 'H';
							break;
						case 'B':
							tmp += 'V';
							break;
						case 'V':
							tmp += 'B';
							break;
						default:
							cerr << "The character (" << c << ") of the PAM is not part of the IUPAC nucleotide nomenclature" << endl;
							tmp += 'N';
							break;
						}

					targetOnDNA[counter_index] = (Tleaf){pamIndices[i], tmp.substr(pamRNA.length()), tmp.substr(pamRNA.length()).c_str(),
									tmp.substr(0, pamRNA.length()),0};
					
				}
			}
			
			counter_index++;
		}
	}
	targetOnDNA.shrink_to_fit();
	
	end = omp_get_wtime();
	cout << end - start << "\n";
	
	cout << "Sorting:\t";
	start = omp_get_wtime(); // sorting the strings before inserting into the tree
	//__gnu_parallel::sort(targetOnDNA.begin(), targetOnDNA.begin() + counter_index, compareFunc);
	sort(targetOnDNA.begin(), targetOnDNA.begin() + counter_index, compareFunc);

	end = omp_get_wtime();
	cout << end - start << "\n";

	//Create tree
	int group_tst = ceil(counter_index / (double)TARG_IN_GROUP); // if a tree is too big, divide it into group_tst smaller trees
	
	for (int jk = 0; jk < group_tst; jk++)
	{
		tree = new Tnode[TARG_IN_GROUP * (pamlen)];
		gruppo = jk+1;
		int inizio = jk * TARG_IN_GROUP;
		int fine = MIN(((jk + 1) * TARG_IN_GROUP), counter_index);
		cout << "Creating tree " << jk + 1 << " of " << group_tst << endl;
		cout << "Build TST:\t";
		
		start = omp_get_wtime(); // build tst
		nodeUsed = 0;
		insall(inizio, fine - 1);
		end = omp_get_wtime();
		cout << end - start << "\n";
		cout << "Save TST:\t";
		start = omp_get_wtime(); // save tst on .bin file
		saveTST(inizio, fine, jk + 1);
		end = omp_get_wtime();
		cout << end - start << "\n";
		delete[] tree;
	}

	globalend = omp_get_wtime(); // end global time
	cout << "-----------------------"
		 << "\n";
	cout << "Total time:\t" << globalend - globalstart << "\n";
	//cout << "C++ end" << endl;
	return 0;
}

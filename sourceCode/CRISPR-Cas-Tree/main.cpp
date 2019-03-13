#include "pamGenerator.cpp"
#include <iostream>
#include <ostream>
#include <omp.h>
#include <string>
#include <vector>
#include <fstream>
#include <parallel/algorithm>

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
	char *guideDNA;
	int next;
} Tleaf;

#define MAXWORDS 32000000
#define MAXCHARS (MAXWORDS * 30)
char space[MAXCHARS];

Tptr tree = (Tptr)malloc(MAXWORDS * 10 * sizeof(Tnode));
int nodeUsed; // number of nodes

Tleaf targetOnDNA[MAXWORDS]; // array of target on DNA
string pamRNA;							 // input pam RNA
string chrName;							 // input chromosome name
string chrSeq;							 // input chromosome sequence

// function that inserts a guide DNA in TST
void insert(char *s, int i)
{
	int d;									 // distance between characters
	s = s + pamRNA.length(); // insert after pam DNA
	Tptr pp = tree;					 // insert start from root
	while (nodeUsed)
	{ // move through TST
		if ((d = *s - pp->splitchar) == 0)
		{
			if (*++s == 0)
			{ // update leaf
				*(s - 1) = 0;
				targetOnDNA[i].next = pp->eqkid;
				pp->eqkid = (i + 1) * -1;
				return;
			}
			if (!pp->eqkid)
			{ // go to eqkid
				pp->eqkid = nodeUsed;
				break;
			}
			pp = &tree[pp->eqkid];
			*(s - 1) = 0;
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
			*(s - 1) = 0;
			pp->eqkid = (i + 1) * -1;
			return;
		}
		*(s - 1) = 0;
		pp->eqkid = nodeUsed;
	}
}

ofstream fileTree;

char pairNuc[2];
bool flag = true;
unsigned char bitNuc;

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
			break; //0 1 0 1   0 0 0 0
		case 'Y':
			bitNuc = 0xA0;
			break; //1 0 1 0   0 0 0 0
		case 'S':
			bitNuc = 0x60;
			break; //0 1 1 0   0 0 0 0
		case 'W':
			bitNuc = 0x90;
			break; //1 0 0 1   0 0 0 0
		case 'K':
			bitNuc = 0xC0;
			break; //1 1 0 0   0 0 0 0
		case 'M':
			bitNuc = 0x30;
			break; //0 0 1 1   0 0 0 0
		case 'B':
			bitNuc = 0xE0;
			break; //1 1 1 0   0 0 0 0
		case 'D':
			bitNuc = 0xD0;
			break; //1 1 0 1   0 0 0 0
		case 'H':
			bitNuc = 0xB0;
			break; //1 0 1 1   0 0 0 0
		case 'V':
			bitNuc = 0x70;
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
			break; //0 0 0 0   0 1 0 1
		case 'Y':
			bitNuc += 0x0A;
			break; //0 0 0 0   1 0 1 0
		case 'S':
			bitNuc += 0x06;
			break; //0 0 0 0   0 1 1 0
		case 'W':
			bitNuc += 0x09;
			break; //0 0 0 0   1 0 0 1
		case 'K':
			bitNuc += 0x0C;
			break; //0 0 0 0   1 1 0 0
		case 'M':
			bitNuc += 0x03;
			break; //0 0 0 0   0 0 1 1
		case 'B':
			bitNuc += 0x0E;
			break; //0 0 0 0   1 1 1 0
		case 'D':
			bitNuc += 0x0D;
			break; //0 0 0 0   1 1 0 1
		case 'H':
			bitNuc += 0x0B;
			break; //0 0 0 0   1 0 1 1
		case 'V':
			bitNuc += 0x07;
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

// Write to file
void saveTST(int arrayDim, int part)
{
	double start, end;

	fileTree.open(pamRNA + "_" + chrName + "_" + to_string(part) + ".bin", ios::out | ios::binary);

	fileTree.write((char *)&arrayDim, sizeof(int)); // write number of targets

	for (int i = 0; i < arrayDim; i++)
	{																																	 // write array of targets on DNA
		fileTree.write((char *)&targetOnDNA[i].guideIndex, sizeof(int)); // write index of target on DNA

		int k = 0;
		bitNuc = 0;
		int counter = 0;
		do
		{ // write target site PAM
			counter++;
			switch (*(targetOnDNA[i].guideDNA))
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
				cerr << "The PAM cannot contain the character (" << *(targetOnDNA[i].guideDNA) << ")" << endl;
				break;
			}
			k++;
			targetOnDNA[i].guideDNA++;
			if (!*(targetOnDNA[i].guideDNA) || k == 2)
			{
				if (counter % 3 == 0 && counter != 0) //TODO da modificare per le pam pari
					bitNuc <<= 4;
				fileTree.put(bitNuc);
				bitNuc = 0;
				k = 0;
			}

			bitNuc <<= 4;
		} while (*(targetOnDNA[i].guideDNA));

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
	serialize(&tree[0]);														// serialize TST

	fileTree.close();
}

/* TIMING */

//Sorting function
bool compareFunc(Tleaf a, Tleaf b)
{
	return strcmp(a.guideDNA + pamRNA.length(), b.guideDNA + pamRNA.length()) < 0;
}

//Insert string in the tree
void insall(int l, int r)
{
	if (r < l)
		return;
	int m = l + (r - l) / 2;
	insert(targetOnDNA[m].guideDNA, m);
	insall(l, m - 1);
	insall(m + 1, r);
}

/**
 * Function that extracts the sequences from the genome. Used only if the user selects the -PAMstart flag
 **/
void extractTargets(vector<int> &pamIndices, int pamlen, int group_tst)
{
	double start, end;
	for (int jk = 0; jk < group_tst; jk++)
	{
		int inizio = jk * 5000000;
		int fine = MIN(((jk + 1) * 5000000), pamIndices.size());
		char *s = space;
		string target;
		cout << "Creating tree " << jk + 1 << " of " << group_tst << endl;
		cout << "Retrieve seq:\t";
		start = omp_get_wtime();
		int i = 0;
		int counter_index = 0;

		for (i = inizio; i < fine; i++)
		{
			if (pamIndices[i] < 0) //entro se pam negativo
			{
				target = chrSeq.substr((pamIndices[i]) * -1, pamlen + 2); //estraggo target+pam dal cromosoma
				if (target.find('N') != std::string::npos)								//se trovo N nel target, salto il target e abbasso il contatore degli indici del cromosoma
				{
					counter_index--;
				}
				else
				{
					targetOnDNA[counter_index] = (Tleaf){pamIndices[i], s, 0}; //salvo l'indice del target
					//cout << "target  <0: " << target << endl;
					//reverse(target.begin(), target.end()); //reverse per aggiungere nell'albero
					for (char &c : target) //salvo la stringa
						*s++ = c;
				}
			}
			else
			{
				target = chrSeq.substr((pamIndices[i]), pamlen + 2);
				//cout << "target >0 prima: " << target << endl;
				if (target.find('N') != std::string::npos)
				{
					counter_index--;
				}
				else
				{
					targetOnDNA[counter_index] = (Tleaf){pamIndices[i], s, 0};
					//cout << "s: ";
					reverse(target.begin(), target.end()); //reverse per aggiungere nell'albero
					//cout << "target  >0: " << target << endl;
					for (char &c : target)
					{ //complemento dei nucleotidi per pam negativa

						switch (c)
						{
						case 'A':
							*s++ = 'T';
							break;
						case 'T':
							*s++ = 'A';
							break;
						case 'C':
							*s++ = 'G';
							break;
						case 'G':
							*s++ = 'C';
							break;
						case 'R':
							*s++ = 'Y';
							break;
						case 'Y':
							*s++ = 'R';
							break;
						case 'S':
							*s++ = 'S';
							break;
						case 'W':
							*s++ = 'W';
							break;
						case 'M':
							*s++ = 'K';
							break;
						case 'K':
							*s++ = 'M';
							break;
						case 'H':
							*s++ = 'D';
							break;
						case 'D':
							*s++ = 'H';
							break;
						case 'B':
							*s++ = 'V';
							break;
						case 'V':
							*s++ = 'B';
							break;
						default:
							cerr << "The character (" << c << ") of the PAM is not part of the IUPAC nucleotide nomenclature" << endl;
							*s++ = 'N';
							break;
						}
						//cout << c;
					}
					// cout << "s: ";
					// for (int jj = 0; jj < 25; jj++){
					// 	cout << s[jj];
					// }
					// cout << endl;
				}
			}
			*s++ = 0;
			counter_index++;
		}

		end = omp_get_wtime();
		cout << end - start << "\n";

		cout << "Sorting:\t";
		start = omp_get_wtime(); // sorting the strins before inserting into the tree
		__gnu_parallel::sort(targetOnDNA, targetOnDNA + counter_index, compareFunc);
		end = omp_get_wtime();
		cout << end - start << "\n";
		cout << "Build TST:\t";
		start = omp_get_wtime(); // build tst
		nodeUsed = 0;
		insall(0, counter_index - 1);
		end = omp_get_wtime();
		cout << end - start << "\n";
		cout << "Save TST:\t";
		start = omp_get_wtime(); // save tst on .bin file
		saveTST(counter_index, jk + 1);
		end = omp_get_wtime();
		cout << end - start << "\n";
	}
}

int main(int argc, char **argv)
{
	string line;																							// line of fasta file
	double start, end, globalstart, globalend, globalpartial; // start and end time, global start and end time
	vector<string> all_pam;																		// vector of all possible pam RNA
	vector<int> pamIndices;																		// vector of target indices of pam RNA on DNA
	ifstream fasta(argv[1]);																	// input fasta file
	ifstream pamfile(argv[2]);																// input pam.txt
	bool variant = atoi(argv[3]);															// variant genome or not
	bool pam_at_start = false;																// pam is at beginning or end of guide
	// if (pam_at_start_s.compare("True") == 0)
	// 	pam_at_start = true;
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
	__gnu_parallel::transform(chrSeq.begin(), chrSeq.end(), chrSeq.begin(), ::toupper); // parallelized to uppercase
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
	int pamlen = pam.length(); //length of the total PAM: (NNNNNNNNNNNNNNNNNNNNNGG) is 23
	if (!pam_at_start)
	{
		pamRNA = pam.substr(pamlen - pamlimit, pamlimit);
	}
	else
	{
		pamRNA = pam.substr(0, pamlimit); // if pam_at_start is set, then PAM = TTTNNNNNNNNNNNNNNNNNNNNN 4, i select the first 4 chars
	}

	all_pam = generatePam(pamRNA, variant); // generate a vector of each possible input pam

	// make list of all possible pam RNA
	string list[all_pam.size()];
	//cout << "allpam size: " << all_pam.size() << ", pam: " << pamRNA << endl;
	for (int j = 0; j < all_pam.size(); j++)
	{
		list[j] = all_pam[j];
		//cout << "pams " << list[j] << endl;
	}

	// ------------------- SEARCH PAM IN THE CHROMOSOME -------------------
	cout << "Search PAM:\t";
	start = omp_get_wtime();
	searchWords(pamIndices, list, all_pam.size(), chrSeq, pamlen, pamlimit, pam_at_start);
	all_pam.clear();
	end = omp_get_wtime();
	cout << end - start << "\n";

	// ------------------------ CREATE THE TST ----------------------------

	// char *t;
	// t = (char *)malloc(80000000 * sizeof(char));
	// free(t);

	double pam_double = pamIndices.size();
	int candidate_targets_per_tst = 10000000;
	int group_tst = ceil(pam_double / candidate_targets_per_tst); // if a tree is too big, divide it into group_tst smaller trees

	if (pam_at_start)
	{
		extractTargets(pamIndices, pamlen, group_tst);
	}
	else
	{ //The user didn't set the -PAMstart flag
		for (int jk = 0; jk < group_tst; jk++)
		{
			int inizio = jk * candidate_targets_per_tst;
			int fine = MIN(((jk + 1) * candidate_targets_per_tst), pamIndices.size());
			char *s = space;
			string target;
			cout << "Creating tree " << jk + 1 << " of " << group_tst << endl;
			cout << "Retrieve seq:\t";
			start = omp_get_wtime();
			int i = 0;
			int counter_index = 0;

			for (i = inizio; i < fine; i++)
			{
				if (pamIndices[i] > 0) //entro se pam positivo
				{
					target = chrSeq.substr(pamIndices[i], pamlen + 2); //estraggo target+pam dal cromosoma
					if (target.find('N') != std::string::npos)				 //se trovo N nel target, salto il target e abbasso il contatore degli indici del cromosoma
					{
						counter_index--;
					}
					else
					{
						targetOnDNA[counter_index] = (Tleaf){pamIndices[i], s, 0}; //salvo l'indice del target

						reverse(target.begin(), target.end()); //reverse per aggiungere nell'albero
						for (char &c : target)								 //salvo la stringa
							*s++ = c;
					}
				}
				else
				{
					target = chrSeq.substr((pamIndices[i]) * -1, pamlen + 2);

					if (target.find('N') != std::string::npos)
					{
						counter_index--;
					}
					else
					{
						targetOnDNA[counter_index] = (Tleaf){pamIndices[i], s, 0};

						for (char &c : target) //complemento dei nucleotidi per pam negativa
							switch (c)
							{
							case 'A':
								*s++ = 'T';
								break;
							case 'T':
								*s++ = 'A';
								break;
							case 'C':
								*s++ = 'G';
								break;
							case 'G':
								*s++ = 'C';
								break;
							case 'R':
								*s++ = 'Y';
								break;
							case 'Y':
								*s++ = 'R';
								break;
							case 'S':
								*s++ = 'S';
								break;
							case 'W':
								*s++ = 'W';
								break;
							case 'M':
								*s++ = 'K';
								break;
							case 'K':
								*s++ = 'M';
								break;
							case 'H':
								*s++ = 'D';
								break;
							case 'D':
								*s++ = 'H';
								break;
							case 'B':
								*s++ = 'V';
								break;
							case 'V':
								*s++ = 'B';
								break;
							default:
								cerr << "The character (" << c << ") of the PAM is not part of the IUPAC nucleotide nomenclature" << endl;
								*s++ = 'N';
								break;
							}
					}
				}
				*s++ = 0;
				counter_index++;
			}

			end = omp_get_wtime();
			cout << end - start << "\n";

			cout << "Sorting:\t";
			start = omp_get_wtime(); // sorting the strins before inserting into the tree
			__gnu_parallel::sort(targetOnDNA, targetOnDNA + counter_index, compareFunc);
			end = omp_get_wtime();
			cout << end - start << "\n";
			cout << "Build TST:\t";
			start = omp_get_wtime(); // build tst
			nodeUsed = 0;
			insall(0, counter_index - 1);
			end = omp_get_wtime();
			cout << end - start << "\n";
			cout << "Save TST:\t";
			start = omp_get_wtime(); // save tst on .bin file
			saveTST(counter_index, jk + 1);
			end = omp_get_wtime();
			cout << end - start << "\n";
		}
	}
	globalend = omp_get_wtime(); // end global time
	cout << "-----------------------"
			 << "\n";
	cout << "Total time:\t" << globalend - globalstart << "\n";
	//cout << "C++ end" << endl;
	return 0;
}
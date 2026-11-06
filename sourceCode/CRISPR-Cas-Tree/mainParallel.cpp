#include "pamGenerator.cpp"
#include "include/indexFormat.h"
#include <iostream>
#include <ostream>
#include <omp.h>
#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <parallel/algorithm>
//#include <algorithm>
//#include <algorithm>
#include <unistd.h>
#include <thread>

/***************************************
 * -18/02/2020
 * 	-Patched PAM at beginning
 * 	-Current Version: YES __gnu_parallel::
 * Note: all databases created before 25/10/2019 are missing the length of the guide in the .bin file
 * 
 */

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace std;

typedef struct tnode *Tptr;
typedef struct tnode
{
	char splitchar;
	int lokid, eqkid, hikid;
} Tnode;

typedef struct tleaf
{
	int guideIndex; // window start ws (signed; see PamSite in pamGenerator.cpp)
	//INSERIRE BITSET CON bitset<#sample> CHE CONTIENE INFO SAMPLES
	string guideDNA;
	string pamDNA;
	int next;
	unsigned char lflags; // CZI_LEAF_PAM_LEFT | CZI_LEAF_FLANK_OFF | pad_valid (HAS_NEXT is added at save)
} Tleaf;

#define MAXWORDS 32000000
#define MAXCHARS (MAXWORDS * 30)

const int TARG_IN_GROUP = 5000000; //number of strings for each TST; Change to have smaller or bigger TSTs

Tptr tree;
int nodeUsed; // number of nodes

vector<Tleaf> targetOnDNA; // array of target on DNA
string pamRNA;			   // input pam RNA
string chrName;			   // input chromosome name
string chrSeq;			   // input chromosome sequence
int gruppo;
bool pam_at_start;
int len_guide_used;

int max_bulges_g;	   // max_bulges, the pad width
bool pam_unconstrained_g; // the build PAM is all N

//the window of a site, pam_len + max_bulges bases from ws, with every position off the
//contig written as 'N'. Off-contig and genomic N are treated alike from here on: in
//the flank both are the 0000 nibble, which matches no PAM letter; in the pad both end
//the bases an alignment may consume; in the guide both discard the site.
string siteWindow(int ws, int width)
{
	long n = (long)chrSeq.size();
	if (ws >= 0 && (long)ws + width <= n)
		return chrSeq.substr(ws, width);
	string w(width, 'N');
	for (int i = 0; i < width; i++)
	{
		long p = (long)ws + i;
		if (p >= 0 && p < n)
			w[i] = chrSeq[p];
	}
	return w;
}

//complement, base by base, of a window (IUPAC codes to their complements, N to N)
string complementSeq(const string &target)
{
	string tmp;
	tmp.reserve(target.size());
	for (char c : target) //complemento dei nucleotidi per pam negativa
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
		case 'N':
			tmp += 'N';
			break;
		case 'Q': //all four bases (include/indexFormat.h): its own complement
			tmp += 'Q';
			break;
		default:
			cerr << "The character (" << c << ") of the PAM is not part of the IUPAC nucleotide nomenclature" << endl;
			tmp += 'N';
			break;
		}
	return tmp;
}

//Applies the format-2 rules to a leaf whose guideDNA (the TST string: the guide read
//from the PAM outward, then the pad) and pamDNA are set. Returns false when the site is
//to be discarded: a hard 'N' (or a base off the contig) in the GUIDE region, which the
//tree has no encoding for. In the pad, the bases before the first 'N' or off-contig
//position are the ones an alignment may consume (pad_valid); the rest are stored as
//'A' -- an alignment that reaches them is rejected at search time whatever they hold,
//so any symbol of the tree's alphabet serves. Until 2.9.1 an 'N' in the pad discarded
//the whole site, though a 0-bulge alignment needs none of it (IDX-13). An 'N' in the
//flank is kept as the 0000 nibble, as in 2.9.1. flank_off: part of the flank lies off
//the contig (only an all-N PAM stores such a site).
bool finishLeaf(Tleaf &leaf, bool left, bool flank_off)
{
	size_t guide_len = leaf.guideDNA.size() - max_bulges_g;
	if (leaf.guideDNA.find('N') < guide_len)
		return false;
	size_t first_n = leaf.guideDNA.find('N', guide_len);
	int pad_valid = first_n == string::npos ? max_bulges_g : (int)(first_n - guide_len);
	for (size_t i = guide_len + pad_valid; i < leaf.guideDNA.size(); i++)
		leaf.guideDNA[i] = 'A';
	leaf.lflags = (left ? CZI_LEAF_PAM_LEFT : 0) | (flank_off ? CZI_LEAF_FLANK_OFF : 0) | (unsigned char)pad_valid;
	return true;
}

// function that inserts a guide DNA in TST
void insert(string st, int i, int i2)
{
	int d; // distance between characters
	const char *s = st.c_str();
	Tptr pp = tree; // insert start from root
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
long num_quad_targets = 0; // stored targets (guide, pad or flank) holding Q
int num_quad_bins = 0;	   // .bin files flagged CZI_FLAG_FOUR_ALLELE

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
		case 'Q': //all four bases, in a root, lo or hi slot (include/indexFormat.h)
			bitNuc = CZI_NIBBLE_Q << 4;
			break; //1 1 1 1   0 0 0 0
		case 'q': //all four bases, in an eq slot, where 1111 ends a string
			bitNuc = CZI_NIBBLE_Q_EQ << 4;
			break; //0 0 0 0   0 0 0 0
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
		case 'Q': //all four bases, in a root, lo or hi slot
			bitNuc += CZI_NIBBLE_Q;
			break; //0 0 0 0   1 1 1 1
		case 'q': //all four bases, in an eq slot
			bitNuc += CZI_NIBBLE_Q_EQ;
			break; //0 0 0 0   0 0 0 0
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

// Serialize TST. Each node's symbol is written in its parent's slot; via_eq says that
// slot is an eq slot, where the four-allele symbol Q is written as 'q' (the 0000
// nibble: 1111 ends a string there, include/indexFormat.h). A tree without Q is
// written byte for byte as by CRISPRitz 2.9.1.
void serialize(Tptr p, bool via_eq)
{
	char sc = (p->splitchar == 'Q' && via_eq) ? 'q' : p->splitchar;
	if (flag)
	{ // write current node and recur for its children
		pairNuc[0] = sc;
		flag = false;
	}
	else
	{
		pairNuc[1] = sc;
		writePair();
		flag = true;
	}
	if (p->lokid > 0)
		serialize(&tree[p->lokid], false); // go to lokid
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
		serialize(&tree[p->hikid], false); // go to hikid
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
		serialize(&tree[p->eqkid], true); // go to eqkid
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

// Whether a leaf's stored strings hold the four-allele symbol Q
bool leafHoldsQ(const Tleaf &leaf)
{
	return leaf.guideDNA.find('Q') != string::npos || leaf.pamDNA.find('Q') != string::npos;
}

// Write to file
void saveTST(int inizio, int fine, int part)
{
	double start, end;
	int arrayDim = fine - inizio;
	// the file is flagged CZI_FLAG_FOUR_ALLELE when one of its leaves holds Q (in the
	// tree string or the flank), so that a reader that cannot decode Q refuses it
	bool holdsQ = false;
	for (int i = inizio; i < fine; i++)
		if (leafHoldsQ(targetOnDNA[i]))
		{
			holdsQ = true;
			num_quad_targets++;
		}
	if (holdsQ)
		num_quad_bins++;
	string fileName = pamRNA + "_" + chrName + "_" + to_string(part) + ".bin";
	fileTree.open(fileName, ios::out | ios::binary);
	if (!fileTree.is_open())
	{
		cerr << "ERROR: cannot create the index file '" << fileName << "'." << endl;
		exit(1);
	}

	// format-2 header (docs/INDEX_FORMAT.md): magic, version, header length, flags,
	// guide length, PAM length, max bulges, then the build PAM region
	int header[6] = {CZI_FORMAT_VERSION, CZI_HEADER_FIXED + (int)pamRNA.size(),
					 (pam_at_start ? CZI_FLAG_PAM_AT_START : 0) | (pam_unconstrained_g ? CZI_FLAG_PAMLESS : 0) |
						 (holdsQ ? CZI_FLAG_FOUR_ALLELE : 0),
					 len_guide_used, (int)pamRNA.size(), max_bulges_g};
	fileTree.write(CZI_MAGIC, 4);
	fileTree.write((char *)header, sizeof(header));
	fileTree.write(pamRNA.data(), pamRNA.size());

	fileTree.write((char *)&arrayDim, sizeof(int));		  // write number of targets
	fileTree.write((char *)&len_guide_used, sizeof(int)); //write len of guide etc etc
	for (int i = inizio; i < fine; i++)
	{																	 // write array of targets on DNA
		fileTree.write((char *)&targetOnDNA[i].guideIndex, sizeof(int)); // write index of target on DNA
		int k = 0;
		bitNuc = 0;
		int counter = 0;
		const char *ppp = targetOnDNA[i].pamDNA.c_str();
		// cout << "pam dna " << targetOnDNA[i].pamDNA << endl;
		do
		{ // write target site PAM
			counter++;
			// cout << "print PAM char to check prima (" << *ppp << ")" << endl;
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
			case 'N':
				bitNuc += 0x0; //'N' is the 0000 nibble: nothing to add
				break;
			case 'Q':
				bitNuc += CZI_NIBBLE_Q; //all four bases: 1111, meets every PAM letter
				break;
			default:
				cerr << "The PAM cannot contain the character (" << *ppp << ")" << endl;
				break;
			}
			k++;
			*ppp++;
			if (!*ppp || k == 2)
			{
				//save two nt in each char then reset the char and k (k keep track of how many nt are already saved in the char [0/1])
				fileTree.put(bitNuc);
				bitNuc = 0;
				k = 0;
			}

			if (k == 1)
				bitNuc <<= 4; //if bitnuc already has one nt written, shift to write the second one
			
		} while (*ppp);

		// the leaf flags (format 2: in place of the '0'/'_' marker), then the next leaf
		// with the same string if there is one
		unsigned char lflags = targetOnDNA[i].lflags | (targetOnDNA[i].next ? CZI_LEAF_HAS_NEXT : 0);
		fileTree.put((char)lflags);
		if (targetOnDNA[i].next)
			fileTree.write((char *)&targetOnDNA[i].next, sizeof(int));
	}
	fileTree.write((char *)&nodeUsed, sizeof(int)); // write number of nodes
	serialize(&tree[0], false);						// serialize TST (the root's slot is no eq slot)
	fileTree.write(CZI_MAGIC, 4);					// end-of-index marker

	// every write above is checked here, at once: a full disk used to leave a cut-short
	// .bin behind an exit 0, which every later search then failed on
	fileTree.close();
	if (fileTree.fail())
	{
		cerr << "ERROR: writing the index file '" << fileName << "' failed (disk full?)." << endl;
		exit(1);
	}
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
	insert(targetOnDNA[m].guideDNA, m, (m - (gruppo - 1) * TARG_IN_GROUP));
	insall(l, m - 1);
	insall(m + 1, r);
}

int main(int argc, char **argv)
{
	//cout << "Versione di test per pam all'inizio -> DONE" << endl;
	string line;											  // line of fasta file
	double start, end, globalstart, globalend, globalpartial; // start and end time, global start and end time
	// vector<string> all_pam;									  // vector of all possible pam RNA
	vector<PamSite> pamIndices;								  // the sites of the contig (window start + shape)
	ifstream fasta(argv[1]);								  // input fasta file
	ifstream pamfile(argv[2]);								  // input pam.txt
	int par_thr = stoi(argv[3]);							  //thread to use in parallel algorithm
	int max_bulges = stoi(argv[4]);							  //max allowed bulges
	int max_mismatches_pam = 0; 							  //future user input for PAM mismatches, currently set to zero
	pam_at_start = false;
	globalstart = omp_get_wtime(); // start global time

	//set parallel threads for gnuparallel algorithms
	omp_set_num_threads(par_thr);

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
	len_guide_used = pamlen - pamlimit;
	
	if (!pam_at_start)
	{
		pamRNA = pam.substr(pamlen - pamlimit, pamlimit);
	}
	else
	{
		pamRNA = pam.substr(0, pamlimit); // if pam_at_start is set, then PAM = TTTNNNNNNNNNNNNNNNNNNNNN -4, i select the first 4 chars
	}

	// a format-2 leaf records how many pad bases are real in 5 bits
	if (max_bulges < 0 || max_bulges > CZI_MAX_BULGES)
	{
		cerr << "ERROR: max bulges must be 0 to " << CZI_MAX_BULGES << " (got " << max_bulges << ")." << endl;
		return 1;
	}
	max_bulges_g = max_bulges;
	pam_unconstrained_g = pamUnconstrained(pamRNA);

	// ------------------- SEARCH PAM IN THE CHROMOSOME -------------------
	start = omp_get_wtime();
	//search pam sequences on the target genome, necessary to extract all putative targets
	pamIndices = searchPAMonGenome(pamRNA, pamlen, chrSeq, pamlimit, pam_at_start, max_bulges, max_mismatches_pam);
	end = omp_get_wtime();
	cout << "Search PAM completed in " << end - start << "\n";

	// ------------------------ CREATE THE TST ----------------------------

	int discarded = 0;
	cout << "Retrieve seq:\t";
	start = omp_get_wtime();
	int i = 0;
	int counter_index = 0;
	targetOnDNA.resize(pamIndices.size());
	const int window = pamlen + max_bulges;
	const int pamWidth = (int)pamRNA.length();
	for (i = 0; i < pamIndices.size(); i++)
	{
		const PamSite &site = pamIndices[i];
		string target = siteWindow(site.ws, window); //pad + guide + PAM, off-contig bases as 'N'
		Tleaf &leaf = targetOnDNA[counter_index];
		leaf.guideIndex = site.ws;
		leaf.next = 0;
		// part of the flank off the contig: left of 0 for a PAM-left window, past the
		// end for a PAM-right one (whose flank follows the pad and the guide)
		bool flank_off = site.left ? site.ws < 0
								   : (long)site.ws + max_bulges + len_guide_used + pamWidth > (long)chrSeq.size();

		// orient the window as 2.9.1 did, into the TST string (guideDNA: the guide from
		// the PAM outward, then the pad) and the stored flank (pamDNA)
		if (pam_at_start && site.left) // 5' PAM, plus strand: [PAM|guide|pad]
		{
			string tmp_pam_str = target.substr(0, pamWidth);
			reverse(tmp_pam_str.begin(), tmp_pam_str.end());
			leaf.guideDNA = target.substr(pamWidth, pamlen - pamWidth + max_bulges);
			leaf.pamDNA = tmp_pam_str;
		}
		else if (pam_at_start) // 5' PAM, minus strand: [pad|guide rc|PAM rc]
		{
			string tmp = complementSeq(target);
			reverse(tmp.begin(), tmp.end());
			string tmp_pam_str = tmp.substr(0, pamWidth);
			reverse(tmp_pam_str.begin(), tmp_pam_str.end());
			leaf.guideDNA = tmp.substr(pamWidth, pamlen - pamWidth + max_bulges);
			leaf.pamDNA = tmp_pam_str;
		}
		else if (!site.left) // 3' PAM, plus strand: [pad|guide|PAM]
		{
			reverse(target.begin(), target.end()); //reverse per aggiungere nell'albero
			leaf.guideDNA = target.substr(pamWidth);
			leaf.pamDNA = target.substr(0, pamWidth);
		}
		else // 3' PAM, minus strand: [PAM rc|guide rc|pad]
		{
			string tmp = complementSeq(target);
			leaf.guideDNA = tmp.substr(pamWidth);
			leaf.pamDNA = tmp.substr(0, pamWidth);
		}

		if (finishLeaf(leaf, site.left, flank_off))
			counter_index++;
		else
			discarded++;
	}
	targetOnDNA.shrink_to_fit();

	end = omp_get_wtime();
	cout << end - start << "\n";

	cout << "Sorting:\t";
	start = omp_get_wtime(); // sorting the strings before inserting into the tree

	//parallel sort of sequences extracted from genome
	__gnu_parallel::sort(targetOnDNA.begin(), targetOnDNA.begin() + counter_index, compareFunc);

	end = omp_get_wtime();
	cout << end - start << "\n";

	//Create tree
	int group_tst = ceil(counter_index / (double)TARG_IN_GROUP); // if a tree is too big, divide it into group_tst smaller trees

	for (int jk = 0; jk < group_tst; jk++)
	{
		tree = new Tnode[TARG_IN_GROUP * (pamlen)];
		gruppo = jk + 1;
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

	if (num_quad_bins)
		cout << "FOUR-ALLELE (Q): " << num_quad_targets << " target(s) hold Q, the symbol of a position "
			 << "whose reference base and SNP ALTs cover all four bases; " << num_quad_bins
			 << " index file(s) flagged for it (docs/INDEX_FORMAT.md)\n";
	globalend = omp_get_wtime(); // end global time
	cout << "-----------------------"
		 << "\n";
	cout << "Total time:\t" << globalend - globalstart << "\n";
	//cout << "C++ end" << endl;
	return 0;
}

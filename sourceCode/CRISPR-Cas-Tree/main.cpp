#include "pamGenerator.cpp"
#include <iostream>
#include <ostream>
#include <omp.h>
#include <string>
#include <vector>
#include <fstream>
#include <parallel/algorithm>

using namespace std;

typedef struct tnode * Tptr;
typedef struct tnode {
    char splitchar;
    int lokid, eqkid, hikid;
} Tnode;

typedef struct tleaf {
    int guideIndex;
	char * guideDNA;
	int next;
} Tleaf;


#define MAXWORDS 25000000
#define MAXCHARS (MAXWORDS*30)
char space[MAXCHARS];

Tptr tree = (Tptr) malloc(MAXWORDS * 10 * sizeof(Tnode));
int nodeUsed;					// number of nodes


Tleaf targetOnDNA[MAXWORDS];	// array of target on DNA
string pamRNA;					// input pam RNA
string chrName;					// input chromosome name
string chrSeq;					// input chromosome sequence




// function that inserts a guide DNA in TST
void insert(char *s, int i){
	int d;								// distance between characters
	s = s + pamRNA.length();			// insert after pam DNA
	Tptr pp = tree;						// insert start from root
	while (nodeUsed) {					// move through TST
        if ((d = *s - pp->splitchar) == 0) {
			if (*++s == 0) {			// update leaf
				*(s-1) = 0;
				targetOnDNA[i].next = pp->eqkid;
				pp->eqkid = (i+1) * -1;
				return;
			}
			if (!pp->eqkid) {			// go to eqkid
				pp->eqkid = nodeUsed;
				break;
			}
			pp = &tree[pp->eqkid];
			*(s-1) = 0;
        } else if (d < 0) {				// go to lowkid
			if (!pp->lokid) {
				pp->lokid = nodeUsed;
				break;
			}
			pp = &tree[pp->lokid];
		} else {
			if (!pp->hikid) {			// go to hikid
				pp->hikid = nodeUsed;
				break;
			}
			pp = &tree[pp->hikid];
		}
    }
    for (;;) {							// insert new node in TST
		pp = &tree[nodeUsed];
        pp->splitchar = *s;
		pp->lokid = pp->eqkid = pp->hikid = 0;
		nodeUsed++;
        if (*++s == 0) {				// add leaf information
			*(s-1) = 0;
			pp->eqkid = (i+1) * -1;
            return;
        }
		*(s-1) = 0;
		pp->eqkid = nodeUsed;
    }
}



ofstream fileTree;


char pairNuc[2];
bool flag = true;
unsigned char bitNuc;

void writePair() {
	
	if (pairNuc[0] == '_')
		bitNuc = 0x30;
	else {
		switch (pairNuc[0]) {
			case 'A': bitNuc = 0x10; break;
			case 'C': bitNuc = 0x20; break;
			case 'G': bitNuc = 0x40; break;
			case 'T': bitNuc = 0x80; break;
			case 'N': bitNuc = 0xF0; break;
			default : bitNuc = 0x0; break;
		}
		switch (pairNuc[1]) {
			case 'A': bitNuc += 0x1; break;
			case 'C': bitNuc += 0x2; break;
			case 'G': bitNuc += 0x4; break;
			case 'T': bitNuc += 0x8; break;
			case 'N': bitNuc += 0xF; break;
			case '_': bitNuc += 0x3; break;
			default : bitNuc += 0x0; break;
		}
	}
	fileTree.put(bitNuc);	
}


// Serialize TST
void serialize(Tptr p) {
	if (flag) {										// write current node and recur for its children
		pairNuc[0] = p->splitchar;
		flag = false;
	} else {
		pairNuc[1] = p->splitchar;
		writePair();
		flag = true;
	}
	if (p->lokid > 0) serialize(&tree[p->lokid]);	// go to lokid
	else {
		if (flag) {
			pairNuc[0] = '0';
			flag = false;
		} else {
			pairNuc[1] = '0';
			writePair();
			flag = true;
		}
	}
	if (p->hikid > 0) serialize(&tree[p->hikid]);	// go to hikid
	else {
		if (flag) {
			pairNuc[0] = '0';
			flag = false;
		} else {
			pairNuc[1] = '0';
			writePair();
			flag = true;
		}
	}
	if (p->eqkid > 0) serialize(&tree[p->eqkid]);	// go to eqkid
	else {
		if (flag)
			pairNuc[0] = '_';
		else
			pairNuc[1] = '_';
		writePair();
		flag = true;
		fileTree.write((char*)&p->eqkid, sizeof(int));	// store array PAM index
	}
	
}

// Write to file
void saveTST(int arrayDim) {
	double start, end;
	//fileTree.open("TST_" + pamRNA + "_" + chrName + ".bin", ios::out | ios::binary);
	fileTree.open(pamRNA + "_" + chrName + ".bin", ios::out | ios::binary);
	
	fileTree.write((char*)&arrayDim, sizeof(int)); 						// write number of targets

	for (int i=0; i<arrayDim; i++) {									// write array of targets on DNA					
		fileTree.write((char*)&targetOnDNA[i].guideIndex, sizeof(int));	// write index of target on DNA

		int k = 0;
		bitNuc = 0;
		do {															// write target site PAM 
			bitNuc <<= 2;
			switch (*(targetOnDNA[i].guideDNA)) {
				case 'A': bitNuc += 0x0; break;
				case 'C': bitNuc += 0x1; break;
				case 'G': bitNuc += 0x2; break;
				case 'T': bitNuc += 0x3; break;
			}
			k++; targetOnDNA[i].guideDNA++;
			if (!*(targetOnDNA[i].guideDNA) || k == 4) {
				fileTree.put(bitNuc);
				bitNuc = 0;
				k = 0;
			}
		} while (*(targetOnDNA[i].guideDNA));
		
		if (targetOnDNA[i].next) {										// write index of next
			fileTree.put('_');
			fileTree.write((char *) &targetOnDNA[i].next, sizeof(int));
		} else {
			fileTree.put('0');
		}
	}
	
	fileTree.write((char*)&nodeUsed, sizeof(int));						// write number of nodes
	serialize(&tree[0]);												// serialize TST		

	fileTree.close();

}



/* TIMING */

bool compareFunc(Tleaf a, Tleaf b) {
    return strcmp(a.guideDNA + pamRNA.length(), b.guideDNA + pamRNA.length()) < 0;
}

void insall(int l, int r) {
	if (r < l) return;
	int m = l + (r - l) / 2;
	insert(targetOnDNA[m].guideDNA, m);
	insall(l, m-1);
	insall(m+1, r);
}


int main( int argc, char **argv) {
	string line;													// line of fasta file
	double start, end, globalstart, globalend, globalpartial;		// start and end time, global start and end time
	vector<string> all_pam;											// vector of all possible pam RNA
	vector<int> pamIndices;											// vector of target indices of pam RNA on DNA
	ifstream fasta(argv[1]);										// input fasta file
	pamRNA = argv[2];												// input PAM
	
	globalstart = omp_get_wtime();									// start global time

// ----------------------- READ INPUT FASTA ---------------------------
	// read chromosome name
	getline(fasta,chrName);
	chrName = chrName.substr(1, chrName.length()-1);
	
	cout << "Load chr:\t"; start = omp_get_wtime();
	while (getline(fasta, line).good())														// read chromosome sequence
		chrSeq+=line;
	__gnu_parallel::transform(chrSeq.begin(), chrSeq.end(), chrSeq.begin(), ::toupper);		// parallelized to uppercase 
	end = omp_get_wtime(); cout << end-start << "\n";

// ------------------- GENERATE ALL POSSIBLE PAM -------------------
	transform(pamRNA.begin(), pamRNA.end(), pamRNA.begin(), ::toupper);	// uppercase of the pam
	all_pam = generatePam(pamRNA);										// generate a vector of each possible input pam
	
	// make list of all possible pam RNA
	string list[all_pam.size()];
	for (int j = 0; j < all_pam.size(); j++) {
		list[j] = all_pam[j];
	}
	
// ------------------- SEARCH PAM IN THE CHROMOSOME -------------------	
	cout << "Search PAM:\t"; start = omp_get_wtime();
	searchWords(pamIndices, list, all_pam.size(), chrSeq);
	all_pam.clear();
	end = omp_get_wtime(); cout << end-start << "\n";

// ------------------------ CREATE THE TST ----------------------------
	char *s = space;
	string target;

	char * t; t = (char *) malloc(80000000*sizeof(char)); free(t);
	
	cout << "Retrieve seq:\t"; start = omp_get_wtime();
	for (int i = 0; i < pamIndices.size(); i++) {  
		if (pamIndices[i] > 0) {
			targetOnDNA[i] = (Tleaf) {pamIndices[i]-1, s, 0};
			target = chrSeq.substr(pamIndices[i]-1, 22 + pamRNA.length());
			reverse(target.begin(), target.end());
			for(char& c : target)
				*s++ = c;
		} else {
			targetOnDNA[i] = (Tleaf) {pamIndices[i]+1, s, 0};
			target = chrSeq.substr((pamIndices[i] + 1)  * -1, 22 + pamRNA.length());
			for(char& c : target)
				switch(c) {
					case 'A': *s++ = 'T'; break;
					case 'T': *s++ = 'A'; break;
					case 'C': *s++ = 'G'; break;
					case 'G': *s++ = 'C'; break;
					default: *s++ = 'N'; break;
				}
		}
		*s++ = 0;
	}
	end = omp_get_wtime(); cout << end-start << "\n";
	
	// parallelized sort by lexicographical order array of target on DNA
	cout << "Sorting:\t"; start = omp_get_wtime();
	__gnu_parallel::sort(targetOnDNA, targetOnDNA + pamIndices.size(), compareFunc);		
	end = omp_get_wtime(); cout << end-start << "\n";
	
	// build TST starting from root
	cout << "Build TST:\t"; start = omp_get_wtime();
	nodeUsed = 0; insall(0, pamIndices.size() - 1);
	end = omp_get_wtime(); cout << end-start << "\n";
	
	// save TST to file
	cout << "Save TST:\t"; start = omp_get_wtime();
	saveTST(pamIndices.size());
	end = omp_get_wtime(); cout << end-start << "\n";
	
	globalend = omp_get_wtime();					// end global time
	cout << "-----------------------" << "\n";
	cout << "Total time:\t" << globalend-globalstart << "\n";
	return 0;
	
}
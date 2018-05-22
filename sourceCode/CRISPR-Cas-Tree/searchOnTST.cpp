#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <string.h>
#include <fstream>
#include <algorithm>

using namespace std;

typedef struct tnode * Tptr;
typedef struct tnode {
    char splitchar;
    int lokid, eqkid, hikid = 0;
} Tnode;

typedef struct tleaf {
    int guideIndex;
	char * guideDNA;
	int next;
} Tleaf;



Tptr tree;
int numNodes;				// number of nodes
Tleaf * targetOnDNA;		// array of target on DNA
int numLeaves;
string pamRNA;				// input pam RNA
string chrName;				// input chromosome name

vector<char *> guideRNA;	// input guide RNA
vector<int> inMM;			// input mismatches
vector<int> inDNAbul;		// input RNA bulge
vector<int> inRNAbul;		// input DNA bulge



ifstream fileTree;

char pairNuc[2];
uint8_t flag = 1;
char in;

void readPair() {
	if (in == 0x30) {
		pairNuc[0] = '_';
		return;
	} else {
		unsigned char first =  in & 0xF0;
		switch (first) {
			case 0x10: pairNuc[0] = 'A'; break;
			case 0x20: pairNuc[0] = 'C'; break;
			case 0x40: pairNuc[0] = 'G'; break;
			case 0x80: pairNuc[0] = 'T'; break;
			case 0xF0: pairNuc[0] = 'N'; break;
			default  : pairNuc[0] = '0'; break;
		}
		unsigned char second = in & 0xF;
		switch (second) {
			case 0x1: pairNuc[1] = 'A'; break;
			case 0x2: pairNuc[1] = 'C'; break;
			case 0x4: pairNuc[1] = 'G'; break;
			case 0x8: pairNuc[1] = 'T'; break;
			case 0xF: pairNuc[1] = 'N'; break;
			case 0x3: pairNuc[1] = '_'; break;
			default : pairNuc[1] = '0'; break;
		}
	}
}


int kid;
int pNode;

void deSerialize() {
	int i = pNode;
	tree[i].splitchar = pairNuc[flag];
	
	if (flag) {
		fileTree.get(in);
		readPair();
		flag = 0;
	} else flag++;	
	if (pairNuc[flag] != '0') { tree[i].lokid = ++pNode; deSerialize(); }	// go to lokid
	
	if (flag) {
		fileTree.get(in);
		readPair();
		flag = 0;
	} else flag++;
	if (pairNuc[flag] != '0') { tree[i].hikid = ++pNode; deSerialize(); }	// go to hikid

	if (flag) {
		fileTree.get(in);
		readPair();
		flag = 0;
	} else flag++;
	if (pairNuc[flag] == '_') { 
		flag++; 
		fileTree.read((char *) &(tree[i].eqkid), sizeof(int)); }
	else { tree[i].eqkid = ++pNode; deSerialize(); }						// go to eqkid
}
	
	
// Load from file by using boost
void loadTST(string path) {	
	fileTree.open(path, ios::in | ios::binary);	
	
	fileTree.read((char*)&numLeaves, sizeof(int));							// read number of leaves
	
	targetOnDNA = (Tleaf *) malloc (numLeaves * sizeof(Tleaf));				// initialize array of targets on DNA
	
	for (int i=0; i<numLeaves; i++) {										// fill array of targets on DNA
		fileTree.read((char*) &targetOnDNA[i].guideIndex, sizeof(int));		// read index of target on DNA
	
		targetOnDNA[i].guideDNA = new char[pamRNA.size()+1];				// initialize PAM size
		targetOnDNA[i].guideDNA[pamRNA.size()] = '\0';
		unsigned char mask;
		int k = 0;
		fileTree.get(in);													// read target PAM
		for (int j=pamRNA.size()-1; j>-1; j--) {
			if (k == 4) {
				fileTree.get(in);
				k = 0;
			}
			mask =  in & 0x3;
			in >>= 2;
			switch (mask) {
				case 0x0: targetOnDNA[i].guideDNA[j] = 'A'; break;
				case 0x1: targetOnDNA[i].guideDNA[j] = 'C'; break;
				case 0x2: targetOnDNA[i].guideDNA[j] = 'G'; break;
				case 0x3: targetOnDNA[i].guideDNA[j] = 'T'; break;
			}
			k++;
		}
		
		fileTree.get(in);													// read index of next PAM with same guide
		//origin.get(in);													// read index of next PAM with same guide
		if (in == '0')
			targetOnDNA[i].next = 0;
		else {
			fileTree.read((char *) &targetOnDNA[i].next , sizeof(int));
		}
		
	}
	

	fileTree.read((char*) &numNodes, sizeof(int));							// read number of nodes
	tree = (Tptr) malloc (numNodes * sizeof(Tnode));						// initialize the TST

	pNode = 0;
	fileTree.get(in);
	readPair();
	flag = 0;
	deSerialize();															// deserialize TST


	fileTree.close();
	
}



char inGuide[23];
int gi;
char targetOfGuide[23];
int ti;

vector<string> vecTargetOfGuide, vecInGuide, bulgeType;
vector<char> directions;
vector<int> indices, mismatches, bulgeSize; 
int bulDNA, bulRNA, mm;



void saveIndices (Tptr p, int d, int bD, int bR, int bulType) {
	if (p->lokid > 0) saveIndices(&tree[p->lokid], d, bD, bR, bulType);		// go to lokid
	if (p->hikid > 0) saveIndices(&tree[p->hikid], d, bD, bR, bulType);		// go to hikid
	
	if (p->eqkid < 0) {											// node is a leaf, save index and strand
		string g(inGuide);										// convert guide to string
		reverse(g.begin(), g.end());
		g += pamRNA;											// add pam to guide
		int index = (p->eqkid + 1) * -1;
		do {
			string t(targetOnDNA[index].guideDNA);				// convert pam to string
			t += targetOfGuide;									// add pam to target
			reverse(t.begin(), t.end());
			vecInGuide.emplace_back(g);							// save guide
			vecTargetOfGuide.emplace_back(t);					// save target
			mismatches.emplace_back(mm-d);						// save mismatches
			if (bulType == 0) {									// NO BULGE case
				bulgeType.emplace_back(" X ");
				bulgeSize.emplace_back(0);
			} else {											// BULGE case
				bulgeType.emplace_back(bulType<0 ? "RNA":"DNA");
				bulgeSize.emplace_back(bulType<0 ? bulRNA-bR:bulDNA-bD);
			}
			if (targetOnDNA[index].guideIndex < 0) {			// negative strand
				indices.emplace_back(targetOnDNA[index].guideIndex * -1);
				directions.emplace_back('-');
			} else {											// strand positive
				indices.emplace_back(targetOnDNA[index].guideIndex + 2 - (bulDNA - bD) + (bulRNA - bR));
				directions.emplace_back('+');
			}
			index = (targetOnDNA[index].next + 1) * -1;
		} while (index > -1);
		
	} else if (p->eqkid > 0) saveIndices(&tree[p->eqkid], d, bD, bR, bulType);	//go to eqkid
}

// CONTROLLARE LA QUESTIONE BULGE CONTEMPORANEI IN DNA ED RNA
void nearsearch(Tptr p, char *s, int d, int bD, int bR, bool goToLoHi, int bulType) {
	
	if (p->lokid > 0 && goToLoHi && (d > 0 || bD > 0 || bR > 0 || *s < p->splitchar))	// go to lokid
		nearsearch(&tree[p->lokid], s, d, bD, bR, true, bulType);
	
    if (p->hikid > 0 && goToLoHi && (d > 0 || bD > 0 || bR > 0 || *s > p->splitchar))	// go to hikid
		nearsearch(&tree[p->hikid], s, d, bD, bR, true, bulType);
	
	if (ti == (20 + bulDNA - bD)) {											// save results
		saveIndices (p, d, bD, bR, bulType);
	} else if (p->eqkid > 0) {												// go to eqkid
		inGuide[gi++] = *s;													// save guide character
		if (*s == p->splitchar) {											// MATCH case
			targetOfGuide[ti++] = p->splitchar;								// save target character uppercase
			nearsearch(&tree[p->eqkid], s+1, d, bD, bR, true, bulType);
			targetOfGuide[ti--] = 0;
		} else if (d > 0) {													// MISMATCH case
			targetOfGuide[ti++] = p->splitchar+32;							// save target character lowercase 
			nearsearch(&tree[p->eqkid], s+1, d-1, bD, bR, true, bulType);
			targetOfGuide[ti--] = 0;
		}
		if (bR > 0 && *(s+1) && bulType < 1) {												// BULGE RNA case
			targetOfGuide[ti++] = '-';										// update last target character with '-'
			nearsearch(p, s+1, d, bD, bR-1, false, bulType-1);
			targetOfGuide[ti--] = 0;
		}
		if (bD > 0 && bulType > -1) {														// BULGE DNA case
			inGuide[gi-1] = '-';											// update last guide character with '-'
			targetOfGuide[ti++] = p->splitchar;								// save target character uppercase
			nearsearch(&tree[p->eqkid], s, d, bD-1, bR, true, bulType+1);
			targetOfGuide[ti--] = 0;
		}
		inGuide[gi--] = 0;
		
	} else if (p->eqkid < 1 && *s) {												// current node is a leaf
		inGuide[gi++] = *s;															// save guide character
		targetOfGuide[ti++] = *s == p->splitchar ? p->splitchar:p->splitchar+32;	// save target character
		d = *s == p->splitchar ? d:d-1;												// update distance
		if (d > -1 && ti == (20 + bulDNA - bD))										// save results
			saveIndices (p, d, bD, bR, bulType);
		d = *s == p->splitchar ? d:d+1;
		targetOfGuide[ti--] = 0;
		inGuide[gi--] = 0;
	}
}





int main( int argc, char **argv) {

	double start, end, globalstart, globalend;		// start and end time, global start and end time	
	string tstFile = argv[1];						// file genome TST
	ifstream fileGuide(argv[2], ios::out);			// file Guide
	mm = atoi(argv[3]); 
	bulDNA = atoi(argv[4]); 
	bulRNA = atoi(argv[5]);
	ofstream fileResults("result.txt", std::ios_base::app);
	int st = tstFile.find_last_of("/");
	int en = tstFile.find_last_of(".");
	pamRNA = tstFile.substr(st, en-st); 
	int un = pamRNA.find_first_of("_");
	chrName = pamRNA.substr(un+1, pamRNA.length());	// retrive chrName
	pamRNA = pamRNA.substr(1, un-1);				// retrive PAM

	globalstart = omp_get_wtime();					// start global time
	
	cout << "Load Guides:\t"; start = omp_get_wtime();
	
	string line;
	int numGuide;
	string iguide;
	
	while (getline(fileGuide, line)) {
		transform(line.begin(), line.end(),line.begin(), ::toupper);	// toUpperCase
		en = line.find_first_of("N");									// find Guide
		if (en > 0) {
			iguide = line.substr(0, en);									// retrive Guide
			reverse(iguide.begin(), iguide.end());
			guideRNA.push_back((char *) malloc(21 * sizeof(char)));
			copy(iguide.begin(), iguide.end(), guideRNA[numGuide]);			// save Guide
			guideRNA[numGuide][20] = '\0';
		/*	line = line.substr(en+1+pamRNA.size(), line.length()-1);		// find mm
			en = line.find_first_of(" ");									// retrive mm
			inMM.push_back(stoi(line.substr(0, en)));						// save mm
			line = line.substr(en+1, line.length()-1); 						// find bulge
			en = line.find_first_of(" ");									// retrive bulge
			inDNAbul.push_back(stoi(line.substr(0, en)));					// save DNA bulge
			inRNAbul.push_back(stoi(line.substr(en, line.length()-1)));		// save RNA bulge*/
		} else {
			reverse(line.begin(), line.end());
			guideRNA.push_back((char *) malloc(21 * sizeof(char)));
			copy(line.begin(), line.end(), guideRNA[numGuide]);			// save Guide
			guideRNA[numGuide][20] = '\0';
		}
		numGuide++;
	}
	end = omp_get_wtime(); cout << end-start << "\n";
	

	// read TST from file
	cout << "Load TST:\t"; start = omp_get_wtime();
	loadTST(tstFile);
	end = omp_get_wtime(); cout << end-start << "\n";

	// Search guide and pam in the TST
	cout << "Nearsearch TST:\t"; start = omp_get_wtime();
	for (int i=0; i<numGuide; i++) {
		ti = gi = 0;
		nearsearch(tree, guideRNA[i], mm, bulDNA, bulRNA, true, 0);
	}
	end = omp_get_wtime(); cout << end-start << "\n";
	
	// Print results
//	cout << "matches = " << indices.size() << endl;
//	cout << "#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge Size" << endl;
	for (int i = 0; i < indices.size(); i++) {
		fileResults << bulgeType[i] << "\t" << vecInGuide[i] << "\t" << vecTargetOfGuide[i] << "\t"<< chrName << "\t" << indices[i] << "\t" << directions[i] << "\t" << mismatches[i] << "\t" << bulgeSize[i] << "\n";
	}
	

	
	globalend = omp_get_wtime();					// end global time
	cout << "-----------------------" << "\n";
	cout << "Total time:\t" << globalend-globalstart << "\n";
	return 0;
	
}
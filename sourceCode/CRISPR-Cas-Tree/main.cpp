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


#define MAXWORDS 29000000
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


//check http://cactus.io/resources/toolbox/decimal-binary-octal-hexadecimal-conversion
void writePair() {
	
	if (pairNuc[0] == '_')
		bitNuc = 0xF0;  //inizialmente 0x30, se lo metto al posto di N 0xF0
	else {
		
		switch (pairNuc[0]) {				//essendo 8 bit, riempio i primi 4
			case 'A': bitNuc = 0x10; break;
			case 'C': bitNuc = 0x20; break;
			case 'G': bitNuc = 0x40; break;
			case 'T': bitNuc = 0x80; break;
			//case 'N': bitNuc = 0xF0; break; //can be deleted
			case 'R': bitNuc = 0x50; break;
			case 'Y': bitNuc = 0xA0; break;
			case 'S': bitNuc = 0x60; break;
			case 'W': bitNuc = 0x90; break;
			case 'K': bitNuc = 0xC0; break;
			case 'M': bitNuc = 0x30; break; //attenzione, è uguale a '_', cambio '_' con N
			case 'B': bitNuc = 0xE0; break;
			case 'D': bitNuc = 0xD0; break;
			case 'H': bitNuc = 0xB0; break;
			case 'V': bitNuc = 0x70; break;
			default : //cout << "pairnuc 0: " << pairNuc[0]<<endl; 
			bitNuc = 0x0; break;
		}
		switch (pairNuc[1]) {				// e poi gli altri 4
			case 'A': bitNuc += 0x1; break;
			case 'C': bitNuc += 0x2; break;
			case 'G': bitNuc += 0x4; break;
			case 'T': bitNuc += 0x8; break;
			//case 'N': bitNuc += 0x0F; break;
			case '_': bitNuc += 0x0F; break;  //inizialmente 0x03, se lo metto al posto di N 0x0F
			case 'R': bitNuc += 0x05; break;
			case 'Y': bitNuc += 0x0A; break;
			case 'S': bitNuc += 0x06; break;
			case 'W': bitNuc += 0x09; break;
			case 'K': bitNuc += 0x0C; break;
			case 'M': bitNuc += 0x03; break; //attenzione, è uguale a '_',  cambio '_' con N
			case 'B': bitNuc += 0x0E; break;
			case 'D': bitNuc += 0x0D; break;
			case 'H': bitNuc += 0x0B; break;
			case 'V': bitNuc += 0x07; break;
			default : //cout << "pairnuc 1: " << pairNuc[1]<<endl;
			bitNuc += 0x0; break;
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
		//cout << "savetst: " << i << "/" <<arrayDim<< endl;
		fileTree.write((char*)&targetOnDNA[i].guideIndex, sizeof(int));	// write index of target on DNA
		//cout << targetOnDNA[i].guideDNA << endl;
		int k = 0;
		bitNuc = 0;
		int counter = 0;
		do {															// write target site PAM 

			counter ++;
			//bitNuc <<= 2;
			//bitNuc <<=4; //TODO è 8 bit, quindi devo sisitemare i numeri
			//cout << "*(targetOnDNA[i].guideDNA " << *(targetOnDNA[i].guideDNA) << endl;
			switch (*(targetOnDNA[i].guideDNA)) {
				case 'A': bitNuc += 0x1; break; // was 0x0
				case 'C': bitNuc += 0x2; break; // was 0x1
				case 'G': bitNuc += 0x4; break; // was 0x2
				case 'T': bitNuc += 0x8; break; // was 0x3
				case 'R': bitNuc += 0x05; break;
				case 'Y': bitNuc += 0x0A; break;
				case 'S': bitNuc += 0x06; break;
				case 'W': bitNuc += 0x09; break;
				case 'K': bitNuc += 0x0C; break;
				case 'M': bitNuc += 0x03; break; //attenzione, è uguale a '_',  cambio '_' con N
				case 'B': bitNuc += 0x0E; break;
				case 'D': bitNuc += 0x0D; break;
				case 'H': bitNuc += 0x0B; break;
				case 'V': bitNuc += 0x07; break;
				default : cerr << "The N char of the pam was not converted";break;//bitNuc +=  0x0; break;

			}
			if (i == 0){
				cout << "bitNuc: ";
                        for (int i = 7; i> -1; i--){
                        auto bit = (bitNuc >> i) & 1U;
                        cout << bit;
                        }
                        cout << endl;
			}

			k++; targetOnDNA[i].guideDNA++;
			if (!*(targetOnDNA[i].guideDNA) || k==2) {//k == 4) {
				//cout<< "in if "<< endl;
				if (counter == 3)
					bitNuc <<=4;
				fileTree.put(bitNuc);
				bitNuc = 0;
				k = 0;
			}
			bitNuc <<=4;
		} while (*(targetOnDNA[i].guideDNA));
		//cout << "end do" << endl; 
		if (targetOnDNA[i].next) {										// write index of next
			//cout << "put _"<<endl;
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
	//cout<<"l "<<l<<" r "<<r<<endl;
	if (r < l) return;
	int m = l + (r - l) / 2;
	insert(targetOnDNA[m].guideDNA, m);
	insall(l, m-1);
	insall(m+1, r);
	//cout<<"faccio insall"<<endl;
}


int main( int argc, char **argv) {
	string line;													// line of fasta file
	double start, end, globalstart, globalend, globalpartial;		// start and end time, global start and end time
	vector<string> all_pam;											// vector of all possible pam RNA
	vector<int> pamIndices;											// vector of target indices of pam RNA on DNA
	ifstream fasta(argv[1]);										// input fasta file
	ifstream pamfile(argv[2]);
   //pamRNA = argv[2];												// input PAM
	
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
	getline(pamfile,line);
   transform(line.begin(), line.end(), line.begin(), ::toupper);// uppercase of the pam
   int delimiter = line.find(" ");
   string pam = line.substr(0, delimiter);
   //cout<<"pam "<<pam<<endl;
   int pamlimit = stoi(line.substr(delimiter, line.length() - 1)); //numero che identifica lunghezza effettiva PAM NNNNNNNNNNNNNNNNNNNNNGG (3)
   int pamlen = pam.length(); //lunghezza della pam totale (NNNNNNNNNNNNNNNNNNNNNGG) 
   pamRNA=pam.substr(pamlen-pamlimit,pamlimit);
	all_pam = generatePam(pam.substr(pamlen-pamlimit,pamlimit));// generate a vector of each possible input pam
	
	// make list of all possible pam RNA
	string list[all_pam.size()];
	
	for (int j = 0; j < all_pam.size(); j++) {
		list[j] = all_pam[j];
      //scout << all_pam[j]<<" "<<endl;
      
	}
	
// ------------------- SEARCH PAM IN THE CHROMOSOME -------------------	
	cout << "Search PAM:\t"; start = omp_get_wtime();
	searchWords(pamIndices, list, all_pam.size(), chrSeq,pamlen,pamlimit);
	all_pam.clear();
	end = omp_get_wtime(); cout << end-start << "\n";
	/*
	for (int l = 0; l < pamIndices.size(); l++){
        cout << "indices " << pamIndices[l] << ","<< endl;
    }*/
// ------------------------ CREATE THE TST ----------------------------	
	char *s = space;
	string target;

	char * t; t = (char *) malloc(80000000*sizeof(char)); free(t);
	
	cout << "Retrieve seq:\t"; start = omp_get_wtime();
	
	//cout << "pam indices size: " << pamIndices.size()<< endl;
	for (int i = 0; i < pamIndices.size(); i++) {
		//cout << "retrieve : " << i << "/" <<pamIndices.size()<< endl;  
		if (pamIndices[i] > 0) {
			target = chrSeq.substr(pamIndices[i], pamlen+2);
			//cout<<"target pos "<<target<<endl;
			if (target.find('N') != std::string::npos){
				//cout << "find n +" << endl;
				//auto it = (pamIndices.begin() + i );
				pamIndices.erase(pamIndices.begin() + i);
				i--;
				continue;
			}
			targetOnDNA[i] = (Tleaf) {pamIndices[i], s, 0};
			//cout << "pam ind [i] -1: " <<pamIndices[i]-1 <<endl;
			//cout << "chr " <<chrSeq.substr(pamIndices[i]-1,1) <<endl;
			//cout << "target + " <<target <<endl;

			reverse(target.begin(), target.end());
			for(char& c : target)
				*s++ = c;
			
		} else {
			target = chrSeq.substr((pamIndices[i])  * -1, pamlen+2);
			//cout<<"target neg "<<target<<endl;
			if (target.find('N') != std::string::npos){
				//cout << "find n -" << endl;
				//auto it = (pamIndices.begin() + i );
				//cout << "Prima "<<*it << endl;
				pamIndices.erase(pamIndices.begin() + i);
				//it = (pamIndices.begin() + i);
				//cout <<"Dopo "<< *it << endl;
				i--;
				continue;
			}//else{
				targetOnDNA[i] = (Tleaf) {pamIndices[i], s, 0};
				
				//cout << "pam ind [i] +1: " << (pamIndices[i] + 1)  * -1<< endl;
				//cout << "chr " <<chrSeq.substr((pamIndices[i] + 1)  * -1,1) <<endl;
				//cout << "target - " <<target <<endl;
				for(char& c : target)
					switch(c) {
						case 'A': *s++ = 'T'; break;
						case 'T': *s++ = 'A'; break;
						case 'C': *s++ = 'G'; break;
						case 'G': *s++ = 'C'; break;
						case 'R': *s++ = 'Y' ; break;
						case 'Y': *s++ = 'R' ; break;
						case 'S': *s++ = 'S' ; break;
						case 'W': *s++ = 'W'; break;
						case 'M': *s++ = 'K'; break;
						case 'K': *s++ = 'M' ; break;
						case 'H': *s++ = 'D' ; break;
						case 'D': *s++ = 'H'; break;
						case 'B': *s++ = 'V'; break;
						case 'V': *s++ = 'B'; break;

						
						default: cout << "c: "<<c<< endl;*s++ = 'N'; break;
					}
			//}
		}
		
		*s++ = 0;
		//cout << "i: " << i << endl;
		//cout << "targ on dna size "<<targetOnDNA[i].guideDNA << endl;
	}

	
	end = omp_get_wtime(); cout << end-start << "\n";
	cout << "size;: " << pamIndices.size()<<endl;
	//cout << "dopo size " << endl;
	// parallelized sort by lexicographical order array of target on DNA
	cout << "Sorting:\t"; start = omp_get_wtime();
	__gnu_parallel::sort(targetOnDNA, targetOnDNA + pamIndices.size(), compareFunc);		
	end = omp_get_wtime(); cout << end-start << "\n";
	//cout << "dopo sorting " << endl;
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

	cout << "C++ end";
	return 0;
	
}
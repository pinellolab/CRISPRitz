#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <string.h>
#include <fstream>
#include <algorithm>

#include <typeinfo> //for printing vars type  typeid(var).name()
#include <bitset>
#include "detailedOutput.h"
#include "convert.h"
#include <stdint.h>
#include <dirent.h> //for directory reading files
#include <array>

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
        //vector<bitset<4>> guideDNA_bit;
        bitset<4> *guideDNA_bit;
} Tleaf;

Tptr tree;
int numNodes;       // number of nodes
Tleaf *targetOnDNA; // array of target on DNA
int numLeaves;
string pamRNA;  // input pam RNA
string chrName; // input chromosome name

vector<char *> guideRNA;                // input guide RNA
vector<string> guideRNA_s;              //input guide RNA as strings;
vector<vector<bitset<4>>> guideRNA_bit; //input guide RNA in bitset

vector<int> inMM;     // input mismatches
vector<int> inDNAbul; // input RNA bulge
vector<int> inRNAbul; // input DNA bulge

ifstream fileTree;

char pairNuc[2];
vector<bitset<4>> pairNuc_bit(2);
uint8_t flag = 1;
char in;

void readPair()
{

        if ((in & 0xF0) == 0xF0)
        { //inizialmente 0x30, se lo metto al posto di N 0xF0
                pairNuc[0] = '_';
                //pairNuc_bit[0] = bitset<4> ("0011");
                pairNuc_bit[0] = in >> 4;
                //cout  << "bit in string uno" << pairNuc_bit[0].to_string() << endl;
                //return;
        }
        else
        {
                unsigned char first = in & 0xF0; // in & 1 1 1 1   0 0 0 0 per tenere solo i primi 4 bit validi (?)
                pairNuc_bit[0] = first >> 4;
                //cout  << "bit in string due" << pairNuc_bit[0].to_string() << endl;
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
                unsigned char second = in & 0xF; //per tenere solo gli ultimi 4 bit validi (?)
                pairNuc_bit[1] = second;
                //cout  << "bit in string " << pairNuc_bit[1].to_string() << endl;
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
                        break; //inizialmente 0x03, se lo metto al posto di N 0x0F
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

int kid;
int pNode;

void deSerialize(vector<Tnode> &albero)
{
        int i = pNode;
        //tree[i].splitchar = pairNuc[flag];
        //tree[i].splitchar_bit = pairNuc_bit[flag];
        albero[i].splitchar_bit = pairNuc_bit[flag];
        if (flag)
        {
                fileTree.get(in);
                readPair();
                flag = 0;
        }
        else
                flag++;
        if (pairNuc_bit[flag] != '0')
        {
                albero[i].lokid = ++pNode;
                deSerialize(albero);
        } // go to lokid

        if (flag)
        {
                fileTree.get(in);
                readPair();
                flag = 0;
        }
        else
                flag++;
        if (pairNuc_bit[flag] != '0')
        {
                albero[i].hikid = ++pNode;
                deSerialize(albero);
        } // go to hikid

        if (flag)
        {
                fileTree.get(in);
                readPair();
                flag = 0;
        }
        else
                flag++;
        if (pairNuc_bit[flag] == '_')
        {
                //if (pairNuc_bit[flag].to_ulong() == 15 || pairNuc_bit[flag].to_ulong() == 240)  {
                flag++;
                fileTree.read((char *)&(albero[i].eqkid), sizeof(int));
        }
        else
        {
                albero[i].eqkid = ++pNode;
                deSerialize(albero);
        } // go to eqkid
}

// Load from file by using boost
//uso 0x0 0x1 0x2 0x3 perchè mi servono solo 2 bit per differenziare le 4 lettere della pam
void loadTST(string path, vector<Tnode> &albero)
{

        fileTree.open(path, ios::in | ios::binary);

        fileTree.read((char *)&numLeaves, sizeof(int));           // read number of leaves
        targetOnDNA = (Tleaf *)malloc(numLeaves * sizeof(Tleaf)); // initialize array of targets on DNA

        for (int i = 0; i < numLeaves; i++)
        { // fill array of targets on DNA
                //cout << "Next" << endl;
                fileTree.read((char *)&targetOnDNA[i].guideIndex, sizeof(int)); // read index of target on DNA
                //targetOnDNA[i].guideDNA_bit.resize(3);
                //targetOnDNA[i].guideDNA = new char[pamRNA.size()+1];        // initialize PAM size
                //targetOnDNA[i].guideDNA[pamRNA.size()] = '\0';
                targetOnDNA[i].guideDNA_bit = new bitset<4>[pamRNA.size()]; //3
                unsigned char mask;
                int k = 0;
                fileTree.get(in);
                int countzero = 0; // read target PAM
                /*
                cout << "bit2: ";
                for (int i = 7; i> -1; i--){
                        auto bit2 = (in >> i) & 1U;
                        cout << bit2;
                }
                cout << endl;
                */
                //cout << "pamrna size: " << pamRNA.size() << endl;  //3
                for (int j = pamRNA.size() - 1; j > -1; j--)
                {
                        if (k == 2)
                        {
                                fileTree.get(in);
                                k = 0;
                        }
                        /*
                        mask =  in & 0x3;
                        in >>= 2; //in = in >> 2
                        */
                        mask = in & 0xF0; //era 0xF
                        in <<= 4;         //era >>=
                        /*
                        cout << "bitX: ";
                        for (int i = 7; i> -1; i--){
                        auto bit = (mask >> i) & 1U;
                        cout << bit;
                        }
                        cout << endl;
                        */
                        //cout << convert(mask) << endl;
                        switch (mask)
                        { //era 0x01 etc
                        case 0x0:
                                countzero++;
                                cout << "Zero" << endl;
                                break;
                        case 0x10:
                                //targetOnDNA[i].guideDNA[j] = 'A';
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0001");
                                //cout << "\n QUI A" << endl;
                                break;
                        case 0x20:
                                //targetOnDNA[i].guideDNA[j] = 'C';
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0010");
                                //cout << "\n QUI C" << endl;
                                break;
                        case 0x40:
                                //targetOnDNA[i].guideDNA[j] = 'G';
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0100");
                                //cout << "\n QUI G" << endl;
                                break;
                        case 0x80:
                                //targetOnDNA[i].guideDNA[j] = 'T';
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1000");
                                //cout << "\n QUI T" << endl;
                                break;
                        case 0x50: //R
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0101");
                                //cout << "\n QUI R" << endl;
                                break;
                        case 0xA0: //Y
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1010");
                                //cout << "\n QUI Y" << endl;
                                break;
                        case 0x60: //S
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0110");
                                //cout << "\n QUI S" << endl;
                                break;
                        case 0x90: //W
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1001");
                                //cout << "\n QUI W" << endl;
                                break;
                        case 0xC0: //K
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1100");
                                //cout << "\n QUI K" << endl;
                                break;
                        case 0x30: //M
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0011");
                                //cout << "\n QUI M" << endl;
                                break;
                        case 0xE0: //B
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1110");
                                //cout << "\n QUI B" << endl;
                                break;
                        case 0xD0: //D
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1101");
                                //cout << "\n QUI D" << endl;
                                break;
                        case 0xB0: //H
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("1011");
                                //cout << "\n QUI H" << endl;
                                break;
                        case 0x70: //V
                                targetOnDNA[i].guideDNA_bit[j] = bitset<4>("0111");
                                //cout << "\n QUI V" << endl;
                                break;
                        default:
                                cerr << "Error on reading the PAM, input char is not an IUPAC nucleotide: " << mask << endl;
                                /*
                                        for (int i = 0; i< 8; i++){
                                                auto bit = (mask >> i) & 1U;
                                                cout << "bit"<<i <<": " << bit << endl;
                                        }
                                        */
                                break;
                        }
                        k++;
                }

                fileTree.get(in); // read index of next PAM with same guide

                //origin.get(in);                             // read index of next PAM with same guide
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

        //tree = (Tptr) malloc (numNodes * sizeof(Tnode));           // initialize the TST
        //tree = new Tnode[numNodes];
        albero.resize(numNodes);
        pNode = 0;
        fileTree.get(in); //creo in usando bitset<4> in modo da prendere esattamente 4 bit, e ogni target sono 23 bitset

        readPair();

        flag = 0;
        deSerialize(albero); // deserialize TST

        fileTree.close();
}

//char inGuide[23];
//vector<bitset<4>> inGuide_bit (23);
//int gi;
//char targetOfGuide[23];
//vector<bitset<4>> targetOfGuide_bit (23);
//int ti;

vector<string> vecTargetOfGuide, vecInGuide, bulgeType;
vector<vector<bitset<4>>> vecTargetOfGuide_bit, vecInGuide_bit; //each target is a vector of bitset
vector<char> directions;
vector<int> indices, mismatches, bulgeSize;
int bulDNA, bulRNA, mm;

void saveIndices(vector<bitset<4>> &inGuide_bit, vector<bitset<4>> &targetOfGuide_bit, vector<Tnode> &p, int pos_tree, int d, int bD, int bR, int bulType, int guideI, vector<vector<int>> &profiling, vector<vector<vector<vector<int>>>> &ext_profiling)
{
        if (p[pos_tree].lokid > 0)
                saveIndices(inGuide_bit, targetOfGuide_bit, p, p[pos_tree].lokid, d, bD, bR, bulType, guideI, profiling, ext_profiling); // go to lokid
        if (p[pos_tree].hikid > 0)
                saveIndices(inGuide_bit, targetOfGuide_bit, p, p[pos_tree].hikid, d, bD, bR, bulType, guideI, profiling, ext_profiling); // go to hikid

        if (p[pos_tree].eqkid < 0)
        { // node is a leaf, save index and strand
                vector<bitset<4>> g_bit(20 + pamRNA.size());

                int j = 0;

                for (auto it = inGuide_bit.rbegin() + pamRNA.size(); it != inGuide_bit.rend(); it++)
                {
                        g_bit[j] = *it;
                        j++;
                }

                //CRITICAL
                for (int i = 0; i < pamRNA.size(); i++)
                {
                        g_bit[g_bit.size() - pamRNA.size() + i] = bitset<4>(convertForBitset(pamRNA[i])); //i can have pam with length > 3
                }

                /* #pragma omp critical
                {
                        for (int i = 0; i < 23; i++)
                                cout << g_bit[i] << endl;
                }*/

                //string g (inGuide);                                       // convert guide to string
                //reverse(g.begin(), g.end());
                //g += pamRNA;
                //cout << typeid(g).name()<< endl;
                //cout<<"inGuide "<<g<<endl;;                    // add pam to guide
                int index = (p[pos_tree].eqkid + 1) * -1;
                do
                {
                        //string t(targetOnDNA[index].guideDNA);                // convert pam to string
                        //t += targetOfGuide; //è il target finale ma al rovescio  // add pam to target
                        //reverse(t.begin(), t.end());

                        vector<bitset<4>> t_bit(20 + pamRNA.size());

                        int i = 0;
                        for (auto it = targetOfGuide_bit.rbegin() + pamRNA.size(); it != targetOfGuide_bit.rend(); it++)
                        {
                                t_bit[i] = *it;
                                i++;
                        }
                        //now t_bit has 20 bitset assigned, need to add the PAM in the last >=3 positions
                        i = t_bit.size() - pamRNA.size();
                        for (int k = pamRNA.size() - 1; k > -1; k--)
                        {
                                t_bit[i] = targetOnDNA[index].guideDNA_bit[k];
                                i++;
                        }

                        //vecInGuide.emplace_back(g);                           // save guide
                        //#pragma omp critical    //faccio un vettore privato e poi li unisco alla fine  in una sezione critica
                        {
                                vecInGuide_bit.emplace_back(g_bit);
                                // posso fare il mio risultato già qui perchè ho già i res di una guida
                                //vecTargetOfGuide.emplace_back(t);                        // save target
                                vecTargetOfGuide_bit.emplace_back(t_bit);

                                mismatches.emplace_back(mm - d); // save mismatches
                                if (bulType == 0)
                                { // NO BULGE case
                                        bulgeType.emplace_back(" X ");
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
                                        directions.emplace_back('-');
                                }
                                else
                                { // strand positive
                                        indices.emplace_back(targetOnDNA[index].guideIndex + 2 - (bulDNA - bD) + (bulRNA - bR));
                                        directions.emplace_back('+');
                                }

                                //    cout<<" indices "<<indices[indices.size()-1]<<endl;
                                index = (targetOnDNA[index].next + 1) * -1;
                                //Profiling e profiling ext
                                detailedOutputFast(guideI, g_bit, t_bit, bulType, mm, profiling, ext_profiling, vecInGuide, vecTargetOfGuide);
                        }
                } while (index > -1);
        }
        else if (p[pos_tree].eqkid > 0)
                saveIndices(inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, d, bD, bR, bulType, guideI, profiling, ext_profiling); //go to eqkid
}

// CONTROLLARE LA QUESTIONE BULGE CONTEMPORANEI IN DNA ED RNA
//void nearsearch(Tptr p, char *s, vector<bitset<4>> guide_bit, int pos_in_guide, int d, int bD, int bR, bool goToLoHi, int bulType, int guideI, vector<vector<int>> &profiling, vector<vector<vector<vector<int>>>> &ext_profiling) { // ,mi son passato la i
void nearsearch(int ti, int gi, vector<bitset<4>> &inGuide_bit, vector<bitset<4>> &targetOfGuide_bit, vector<Tnode> &p, int pos_tree, vector<bitset<4>> guide_bit, int pos_in_guide, int d, int bD, int bR, bool goToLoHi, int bulType, int guideI, vector<vector<int>> &profiling, vector<vector<vector<vector<int>>>> &ext_profiling)
{ // ,mi son passato la i

        // uso hash bitset al posto della conversione
        if (p[pos_tree].lokid > 0 && goToLoHi && (d > 0 || bD > 0 || bR > 0 || guide_bit[pos_in_guide].to_ulong() < p[pos_tree].splitchar_bit.to_ulong())) // go to lokid
                nearsearch(ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].lokid, guide_bit, pos_in_guide, d, bD, bR, true, bulType, guideI, profiling, ext_profiling);

        if (p[pos_tree].hikid > 0 && goToLoHi && (d > 0 || bD > 0 || bR > 0 || guide_bit[pos_in_guide].to_ulong() > p[pos_tree].splitchar_bit.to_ulong())) // go to hikid
                nearsearch(ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].hikid, guide_bit, pos_in_guide, d, bD, bR, true, bulType, guideI, profiling, ext_profiling);

        if (ti == (20 + bulDNA - bD))
        {                                                                                                                       // save results
                saveIndices(inGuide_bit, targetOfGuide_bit, p, pos_tree, d, bD, bR, bulType, guideI, profiling, ext_profiling); //salvo la i
        }
        else if (p[pos_tree].eqkid > 0)
        { // go to eqkid

                //inGuide[gi++] = *s;                                        // save guide character
                //BEGIN CRITICAL
                inGuide_bit[gi] = guide_bit[pos_in_guide]; //remember add - 1
                //cout << convert(guide_bit[pos_in_guide]) <<", thr num: " << omp_get_thread_num() << endl;
                //cout << convert(guide_bit[pos_in_guide]) << " " ;
                gi++;
                //END CRITICAL
                //if (*s == p->splitchar) {   // MATCH case
                if ((guide_bit[pos_in_guide] & p[pos_tree].splitchar_bit) != 0)
                { //MATCH CASE

                        //targetOfGuide[ti++] = p->splitchar;                     // save target character uppercase
                        targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit; //remember add -1
                        ti++;
                        pos_in_guide++;
                        nearsearch(ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, guide_bit, pos_in_guide, d, bD, bR, true, bulType, guideI, profiling, ext_profiling);
                        ti--;
                        //targetOfGuide[ti--] = 0;
                }
                else if (d > 0)
                { // MISMATCH case
                        //targetOfGuide[ti++] = p->splitchar+32;                       // save target character lowercase
                        targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit; //remeber add -1
                        ti++;
                        pos_in_guide++;
                        nearsearch(ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, guide_bit, pos_in_guide, d - 1, bD, bR, true, bulType, guideI, profiling, ext_profiling);
                        ti--;
                        //targetOfGuide[ti--] = 0;
                }
                //if (bR > 0 && *(s+1) && bulType < 1) {    // BULGE RNA case
                if (bR > 0 && (pos_in_guide + 1) < 20 && bulType < 1)
                { // BULGE RNA case
                        //targetOfGuide[ti++] = '-';                                       // update last target character with '-'
                        targetOfGuide_bit[ti] = bitset<4>("0000"); // il char '-'
                        ti++;
                        pos_in_guide++;
                        nearsearch(ti, gi, inGuide_bit, targetOfGuide_bit, p, pos_tree, guide_bit, pos_in_guide, d, bD, bR - 1, false, bulType - 1, guideI, profiling, ext_profiling);
                        //targetOfGuide[ti--] = 0;
                        ti--;
                }
                if (bD > 0 && bulType > -1)
                { // BULGE DNA case
                        //inGuide[gi-1] = '-';                                           // update last guide character with '-'
                        inGuide_bit[gi - 1] = bitset<4>("0000"); //il char '-'
                        //targetOfGuide[ti++] = p[pos_tree].splitchar;                             // save target character uppercase
                        targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit;
                        ti++;
                        nearsearch(ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, guide_bit, pos_in_guide, d, bD - 1, bR, true, bulType + 1, guideI, profiling, ext_profiling);
                        //targetOfGuide[ti--] = 0;
                        ti--;
                }
                //inGuide[gi--] = 0;
                gi--;

                //} else if (p->eqkid < 1 && *s) {                                                // current node is a leaf
        }
        else if (p[pos_tree].eqkid < 1 && pos_in_guide < 20)
        { //current node is a leaf

                //inGuide[gi++] = *s;                                                      // save guide character
                //BEGIN CRITICAL
                inGuide_bit[gi] = guide_bit[pos_in_guide]; //rememeber add -1
                gi++;
                //END CRITICAL
                //targetOfGuide[ti++] = *s == p->splitchar ? p->splitchar:p->splitchar+32;        // save target character       d = *s == p->splitchar ? d:d-1;                                                                                         // update distance

                targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit; //remebere add -1
                ti++;

                if (d > -1 && ti == (20 + bulDNA - bD))
                { // save results
                        saveIndices(inGuide_bit, targetOfGuide_bit, p, pos_tree, d, bD, bR, bulType, guideI, profiling, ext_profiling);
                }
                d = (guide_bit[pos_in_guide] & p[pos_tree].splitchar_bit) == 0 ? d : d + 1;
                //targetOfGuide[ti--] = 0;
                ti--;
                gi--;
                //inGuide[gi--] = 0;
        }
}

void resetGlbVar()
{
        double start, end; //for timing
        cout << "Reset time:\t";
        start = omp_get_wtime();

        numNodes = 0; // number of nodes

        for (int i = 0; i < numLeaves; i++)
        {
                delete[] targetOnDNA[i].guideDNA_bit;
        }

        free(targetOnDNA); // array of target on DNA

        //delete tree;
        numLeaves = 0;

        inMM.clear();     // input mismatches
        inDNAbul.clear(); // input RNA bulge
        inRNAbul.clear(); // input DNA bulge
        //pairNuc.clear(); pairNuc.resize(2);
        pairNuc_bit.clear();
        pairNuc_bit.resize(2);
        flag = 1;
        //in = '0';
        kid = 0;
        pNode = 0;
        //char inGuide[23]; // non servirà
        //inGuide_bit.clear(); inGuide_bit.resize(23);
        //char targetOfGuide[23]; // non servirà
        //targetOfGuide_bit.clear(); targetOfGuide_bit.resize(23);

        vecTargetOfGuide.clear();
        vecInGuide.clear();
        bulgeType.clear();

        vecTargetOfGuide_bit.clear(), vecInGuide_bit.clear(); //each target is a vector of bitset

        directions.clear();
        indices.clear();
        mismatches.clear();
        bulgeSize.clear();

        end = omp_get_wtime();
        cout << end - start << "\n";
}

int test_for_omp;

int main(int argc, char **argv)
{
        DIR *d;
        dirent *dir;
        vector<string> fileList;
        string line;
        double start, end, globalstart, globalend; // start and end time, global start and end time
        //string tstFile = argv[1];                               // file genome TST

        ifstream fileGuide(argv[2]); // file Guide
        mm = atoi(argv[3]);
        bulDNA = atoi(argv[4]);
        bulRNA = atoi(argv[5]);
        ifstream pamfile(argv[6]);
        char *genome_dir = argv[1];
        getline(pamfile,line);
        transform(line.begin(), line.end(), line.begin(), ::toupper);// uppercase of the pam
        int delimiter = line.find(" ");
        string pam = line.substr(0, delimiter);
        //cout<<"pam "<<pam<<endl;
        int pamlimit = stoi(line.substr(delimiter, line.length() - 1)); //numero che identifica lunghezza effettiva PAM NNNNNNNNNNNNNNNNNNNNNGG (3)
        int pamlen = pam.length(); //lunghezza della pam totale (NNNNNNNNNNNNNNNNNNNNNGG) 
        pamRNA=pam.substr(pamlen-pamlimit,pamlimit);
        globalstart = omp_get_wtime(); // start global time

        cout << "Load Guides:\t";
        //start = omp_get_wtime();

        int numGuide;
        string iguide;
        int en;
        while (getline(fileGuide, line))
        {
                transform(line.begin(), line.end(), line.begin(), ::toupper); // toUpperCase
                en = line.find_first_of("N");                                 // find Guide
                if (en > 0)
                {
                        
                        iguide = line.substr(0, en); // retrive Guide
                        guideRNA_s.push_back(iguide);
                        reverse(iguide.begin(), iguide.end());
                        guideRNA.push_back((char *)malloc(21 * sizeof(char)));
                        copy(iguide.begin(), iguide.end(), guideRNA[numGuide]); // save Guide
                        guideRNA[numGuide][20] = '\0';
                        /*      line = line.substr(en+1+pamRNA.size(), line.length()-1);     // find mm
                        en = line.find_first_of(" ");                                            // retrive mm
                        inMM.push_back(stoi(line.substr(0, en)));                               // save mm
                        line = line.substr(en+1, line.length()-1);                            // find bulge
                        en = line.find_first_of(" ");                                        // retrive bulge
                        inDNAbul.push_back(stoi(line.substr(0, en)));                       // save DNA bulge
                        inRNAbul.push_back(stoi(line.substr(en, line.length()-1)));         // save RNA bulge*/
                }
                else
                {
                        guideRNA_s.push_back(line);
                        reverse(line.begin(), line.end());
                        guideRNA.push_back((char *)malloc(21 * sizeof(char)));
                        copy(line.begin(), line.end(), guideRNA[numGuide]); // save Guide
                        guideRNA[numGuide][20] = '\0';
                }
                numGuide++;
        }

        //transform loaded guides into bitset //TODO insert into previous for loop
        for (int i = 0; i < guideRNA.size(); i++)
        {
                vector<bitset<4>> tmp(20);
                for (int j = 0; j < 20; j++)
                {
                        tmp[j] = bitset<4>(convert(guideRNA[i][j]));
                }
                guideRNA_bit.push_back(tmp);
        }

        end = omp_get_wtime();
        cout << end - start << "\n";

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

                        fileList.push_back(dir->d_name);
                }
        }
        closedir(d);
        ofstream fileResults("result.txt", std::ios_base::app);
        ofstream fileprofiling("profile_test_fast.xls", std::ios_base::out); //also add specific name //
        ofstream file_ext_profiling("profile_ext_test_fast.xls", std::ios_base::out);
        //Create matrix for profiling
        vector<vector<int>> profiling(22 + mm + 1, vector<int>(numGuide, 0)); //22+mm+1 is number of cols in a profile
        //TODO controllare i numeri
        //Cols labels for profiling output file
        fileprofiling << "GUIDE\t";
        for (int i = 0; i < 20; i++)
        {
                fileprofiling << "BP\t";
        }
        fileprofiling << "\t"
                      << "ONT\tOFFT\tOFFT/MM\t\t";
        for (int i = 0; i <= mm; i++)
        {
                fileprofiling << i << "MM\t";
        }
        fileprofiling << "\n";

        //Create matrix for extended profiling
        vector<vector<vector<vector<int>>>> ext_profiling(numGuide, vector<vector<vector<int>>>(mm + 1, vector<vector<int>>(4, vector<int>(20, 0)))); //TODO thiss
        //1 dim = id of guide, 2 dim = select matrix with an x number of mms in the target, 3 dim = select nucleotide (acgt), 4 dim = select position

        //inGuide_bit.resize(20 + pamRNA.size());
        //targetOfGuide_bit.resize(20 + pamRNA.size());
        //FOR ALL THE CHR
        for (int file = 0; file < fileList.size(); file++)
        {

                cout << "File " << fileList[file] << endl;
                string gen_dir(genome_dir);
                int point = fileList[file].find_last_of(".");
                int under = fileList[file].find_first_of("_");

                chrName = fileList[file].substr(under + 1, point - under - 1);

                //cout << "Pam rna " << pamRNA << ", chrname: " << chrName << endl;

                //int underscore = fileList[file].find_first_of("_");

                //int st = fileList[file].find_last_of("/");
                //int st = 0;

                //en = fileList[file].find_last_of(".");

                //cout << "st: " << st << ", en: " << en << endl;
                //pamRNA = fileList[file].substr(st, en-st);

                //int un = pamRNA.find_first_of("_");
                //chrName = pamRNA.substr(un+1, pamRNA.length()); // retrive chrName
                //pamRNA = pamRNA.substr(1, un-1);                                // retrive PAM

                //cout << "PAM rna " << pamRNA << endl;

                // create tree
                vector<Tnode> albero;

                // read TST from file
                cout << "Load TST:\t";
                start = omp_get_wtime();
                string aa = gen_dir + "/" + fileList[file];
                cout << "carico " << aa << endl;
                loadTST(gen_dir + "/" + fileList[file], albero);

                end = omp_get_wtime();
                cout << end - start << "\n";

                cout << "Albero size: " << albero.size() << endl;

                /*
                for(int i = 0; i < numLeaves; i++){
                        cout << "Pam: ";
                        for (int j = 0; j< 3 ; j++)
                             cout<<   convert(targetOnDNA[i].guideDNA_bit[j]);
                        cout << endl;
                }
                */
                // Search guide and pam in the TST
                cout << "Nearsearch TST:\t";
                start = omp_get_wtime();

                double newguide = numGuide;
                int group_guide = ceil(newguide / 1000);

                for (int jk = 0; jk < group_guide; jk++)
                {

                        int inizio = jk * 1000;
                        int fine = MIN(((jk + 1) * 1000), numGuide);
                        int m_t = omp_get_max_threads();
                        //m_t = 1;
                        int mm2 = mm;
                        int bulDNA2 = bulDNA;
                        int bulRNA2 = bulRNA;
                        //#pragma omp parallel  num_threads(2) firstprivate(mm2,bulDNA2, bulRNA2)
                        {
                                //Tolto targetOnDna
                                vector<bitset<4>> inGuide_bit(23); //mettere fuori e dichiararle private
                                vector<bitset<4>> targetOfGuide_bit(23);
                                //#pragma omp for
                                for (int i = inizio; i < fine; i++)
                                {
                                        int ti = 0;
                                        int gi = 0;
                                        //cout << "\nguide RNA: " << guideRNA[i] << endl;
                                        //#pragma omp single
                                        {
                                                //#pragma omp task
                                                {
                                                        nearsearch(ti, gi, inGuide_bit, targetOfGuide_bit, albero, 0, guideRNA_bit[i], 0, mm2, bulDNA2, bulRNA2, true, 0, i, profiling, ext_profiling); //mi passo la i;
                                                }
                                        }
                                        //CONTROLLO
                                        //inGuide_bit.clear(); inGuide_bit.resize(23);
                                        //targetOfGuide_bit.clear(); targetOfGuide_bit.resize(23);
                                }
                        }
                        cout << "indices size: " << indices.size() << endl;
                        for (int i = 0; i < indices.size(); i++)
                        {
                                fileResults << bulgeType[i] << "\t" << vecInGuide[i] << "\t" << vecTargetOfGuide[i] << "\t" << chrName << "\t" << indices[i] << "\t" << directions[i] << "\t" << mismatches[i] << "\t" << bulgeSize[i] << "\n";
                                //cout << bulgeType[i] << "\t" << vecInGuide[i] << "\t" << vecTargetOfGuide[i] << "\t"<< chrName << "\t" << indices[i] << "\t" << directions[i] << "\t" << mismatches[i] << "\t" << bulgeSize[i] << "\n";
                        }

                        //Test duplicati

                        sort(indices.begin(), indices.end());
                        int dup_count = 0;
                        // for (int i = 0; i < indices.size(); i++)
                        // {
                        //         if (indices[i] == indices[i - 1])
                        //         {
                        //                 cout << "duplicato: " << indices[i] << endl;
                        //                 dup_count++;
                        //         }
                        // }
                        // cout << "Num duplicati: " << dup_count << endl;

                        bulgeType.clear();
                        vecInGuide.clear();
                        vecTargetOfGuide.clear();
                        indices.clear();
                        directions.clear();
                        mismatches.clear();
                        bulgeSize.clear();
                }

                end = omp_get_wtime();
                cout << end - start << "\n";
                globalend = omp_get_wtime(); // end global time
                cout << "-----------------------"
                     << "\n";
                cout << "Total time:\t" << globalend - globalstart << "\n";

                //reset all global variables for next chr

                resetGlbVar();
        }
        for (int i = 0; i < numGuide; i++)
        {
                saveProfileGuide(guideRNA_s[i], i, mm, profiling, ext_profiling, fileprofiling, file_ext_profiling);
        }
        fileprofiling.close();
        file_ext_profiling.close();

        cout << "\nC++ end" << endl;
        return 0;
}
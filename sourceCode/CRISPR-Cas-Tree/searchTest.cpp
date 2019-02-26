#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <algorithm>

#include <typeinfo> //for printing vars type  typeid(var).name()
#include <bitset>
#include "detailedOutput.h"
#include "convert.h"
#include <stdint.h>
#include <dirent.h> //for directory reading files
#include <array>
#include <sys/stat.h>

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

Tptr tree;
//int numNodes;       // number of nodes
//Tleaf *targetOnDNA; // array of target on DNA
//int numLeaves;
string pamRNA;  // input pam RNA

vector<char *> guideRNA;                // input guide RNA
vector<string> guideRNA_s;              //input guide RNA as strings;
vector<vector<bitset<4>>> guideRNA_bit; //input guide RNA in bitset





// vector<bitset<4>> pairNuc_bit(2);
// char pairNuc[2];
// uint8_t flag = 1;
// char in;
int len_guide;

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
                unsigned char first = in & 0xF0;        // keep only 4 left-most bits
                pairNuc_bit[0] = first >> 4;            // shift them to the right to save them in a bitset<4>
                switch (first) {
                        case 0x10: pairNuc[0] = 'A'; break;
                        case 0x20: pairNuc[0] = 'C'; break;
                        case 0x40: pairNuc[0] = 'G'; break;
                        case 0x80: pairNuc[0] = 'T'; break;
                        //case 0xF0: pairNuc[0] = 'N'; break;
                        case 0x50: pairNuc[0] = 'R'; break;
                        case 0xA0: pairNuc[0] = 'Y'; break;
                        case 0x60: pairNuc[0] = 'S'; break;
                        case 0x90: pairNuc[0] = 'W'; break;
                        case 0xC0: pairNuc[0] = 'K'; break;
                        case 0x30: pairNuc[0] = 'M'; break;
                        case 0xE0: pairNuc[0] = 'B'; break;
                        case 0xD0: pairNuc[0] = 'D'; break;
                        case 0xB0: pairNuc[0] = 'H'; break;
                        case 0x70: pairNuc[0] = 'V'; break;
                        default  : pairNuc[0] = '0'; break;
                }
                unsigned char second = in & 0xF;        // keep only 4 right-most bits
                pairNuc_bit[1] = second;
                switch (second) {
                        case 0x1: pairNuc[1] = 'A'; break;
                        case 0x2: pairNuc[1] = 'C'; break;
                        case 0x4: pairNuc[1] = 'G'; break;
                        case 0x8: pairNuc[1] = 'T'; break;
                        //case 0xF: pairNuc[1] = 'N'; break;
                        case 0x0F: pairNuc[1] = '_'; break; //inizialmente 0x03, se lo metto al posto di N 0x0F
                        case 0x05: pairNuc[1] = 'R'; break;
                        case 0x0A: pairNuc[1] = 'Y'; break;
                        case 0x06: pairNuc[1] = 'S'; break;
                        case 0x09: pairNuc[1] = 'W'; break;
                        case 0x0C: pairNuc[1] = 'K'; break;
                        case 0x03: pairNuc[1] = 'M'; break;
                        case 0x0E: pairNuc[1] = 'B'; break;
                        case 0x0D: pairNuc[1] = 'D'; break;
                        case 0x0B: pairNuc[1] = 'H'; break;
                        case 0x07: pairNuc[1] = 'V'; break;
                        default :  pairNuc[1] = '0'; break;
                }
        }
}

int kid;
//int pNode;

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
                deSerialize(albero,fileTree, pairNuc_bit, pairNuc, in, flag, pNode);
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
Tleaf* loadTST(string path, vector<Tnode> &albero, ifstream &fileTree, int &numNodes, int &numLeaves )
{       
        Tleaf * targetOnDNA;
        fileTree.open(path, ios::in | ios::binary);

        fileTree.read((char *)&numLeaves, sizeof(int));           // read number of leaves
        targetOnDNA = (Tleaf *)malloc(numLeaves * sizeof(Tleaf)); // initialize array of targets on DNA
        char in;
        for (int i = 0; i < numLeaves; i++)
        { // fill array of targets on DNA
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
                        
                        /*
                        cout << "bitX: ";
                        for (int i = 7; i> -1; i--){
                        auto bit = (mask >> i) & 1U;
                        cout << bit;
                        }
                        cout << endl;
                        */

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

//char inGuide[23];
//vector<bitset<4>> inGuide_bit (23);
//int gi;
//char targetOfGuide[23];
//vector<bitset<4>> targetOfGuide_bit (23);
//int ti;

// vector<string> vecTargetOfGuide, vecInGuide, bulgeType;
// vector<vector<bitset<4>>> vecTargetOfGuide_bit, vecInGuide_bit; //each target is a vector of bitset
// vector<char> directions;
// vector<int> indices, mismatches, bulgeSize;
vector <string> save_glb;
int bulDNA, bulRNA, mm;

bool create_target;
bool create_profile;

void saveIndices( vector<bitset<4>> &inGuide_bit, vector<bitset<4>> &targetOfGuide_bit, vector<Tnode> &p, int pos_tree, int d, int bD, int bR, int bulType, 
                Tleaf *targetOnDNA, char (&inGuide)[23], char (&targetOfGuide)[23], vector<string> &vecTargetOfGuide, vector<string> &vecInGuide, vector<string> &bulgeType, vector<char> &directions, vector<int> &indices,  vector<int> &mismatches,  vector<int> &bulgeSize,
                int guideI, vector<vector<vector<int>>> &profiling, vector<vector<vector<vector<vector<int>>>>> &ext_profiling, vector<vector<vector<int>>> &pd, vector<vector<vector<int>>> &pdm, vector<vector<vector<int>>> &pr, vector<vector<vector<int>>> &prm)
{//char (&inGuide)[23], char (&targetOfGuide)[23]
        if (p[pos_tree].lokid > 0)
                saveIndices( inGuide_bit, targetOfGuide_bit, p, p[pos_tree].lokid, d, bD, bR, bulType, targetOnDNA, 
                inGuide, targetOfGuide, vecTargetOfGuide, vecInGuide,bulgeType, directions,indices, mismatches, bulgeSize,
                guideI, profiling, ext_profiling, pd, pdm, pr, prm); // go to lokid
        if (p[pos_tree].hikid > 0)
                saveIndices( inGuide_bit, targetOfGuide_bit, p, p[pos_tree].hikid, d, bD, bR, bulType, targetOnDNA, 
                inGuide, targetOfGuide, vecTargetOfGuide, vecInGuide,bulgeType, directions,indices, mismatches, bulgeSize,
                guideI, profiling, ext_profiling, pd, pdm, pr, prm); // go to hikid

        if (p[pos_tree].eqkid < 0)
        { // node is a leaf, save index and strand
                //convert guide to string
                string g(inGuide); // convert guide to string
                //cout << "in guide: " << inGuide << endl;
                reverse(g.begin(), g.end());
                g += pamRNA;
                
                
                vector<bitset<4>> g_bit(len_guide + pamRNA.size() + bulDNA);
                
                int j = 0;
                
                for (j = 0; j < (len_guide+bulDNA - bD); j++){          
                        g_bit[j] = inGuide_bit[j];
                }

                reverse(g_bit.begin(), g_bit.begin() + len_guide + bulDNA - bD);
                
                int pos_to_save = 0;
                for (int i = 0; i < pamRNA.size(); i++)         
                {
                        g_bit[g_bit.size() - pamRNA.size() - bD + i] = bitset<4>(convertCharToBitset(pamRNA[i]));
                        pos_to_save = g_bit.size() - pamRNA.size() -bD + i;
                }
                pos_to_save++;
                g_bit.resize(pos_to_save);

                int index = (p[pos_tree].eqkid + 1) * -1;
                // //Test for dup
                // string guide_correct_dir;
                // for (int i = 0; i < g_bit.size(); i++){
                //         guide_correct_dir += convertBitsetToChar(g_bit[i]);
                // }
                // duplicati << guide_correct_dir << "\t" << index << "\n";
                do
                {
                        // string t(targetOnDNA[index].guideDNA); // convert pam to string
                        // t += targetOfGuide;                    // add pam to target
                        // reverse(t.begin(), t.end());
                        string t (targetOfGuide);
                        //cout << "target: " << targetOfGuide << endl;
                        reverse(t.begin(),t.end());
                        t+= targetOnDNA[index].guideDNA;

                        vector<bitset<4>> t_bit(len_guide + pamRNA.size() + bulDNA);
                        
                        int i = 0;
                        for (i = 0; i < (len_guide+bulDNA - bD) ; i++){
                                t_bit[i] = targetOfGuide_bit[i];
                        }


                        reverse(t_bit.begin(), t_bit.begin() + len_guide + bulDNA - bD);
                      
                        i = len_guide + bulDNA -bD;
                        
                        for (int k = 0; k < pamRNA.size();k++)
                        {
                                t_bit[i] = targetOnDNA[index].guideDNA_bit[k];
                                i++;
                        }
                        
                        t_bit.resize(pos_to_save);
                       

                        //Controllo che tipo di result voglio
                        if(create_target){
                                
                                
                                //vecInGuide_bit.push_back(g_bit);
                                // posso fare il mio risultato già qui perchè ho già i res di una guida
                                //vecTargetOfGuide_bit.emplace_back(t_bit);

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

                                //Profiling e profiling ext
                                if (create_profile){
                                        if (bulType == 0)
                                                detailedOutputFast(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, profiling, ext_profiling, vecInGuide, vecTargetOfGuide);
                                        else if (bulType > 0)
                                                detailedOutputFastBulgeDNA(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, bD, bulDNA, pd, pdm, vecInGuide, vecTargetOfGuide);
                                        else
                                                detailedOutputFastBulgeRNA(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, bD, bulDNA, pr, prm, vecInGuide, vecTargetOfGuide);
                                } else{
                                        // string guide_correct_dir;
                                        // string target_correct_dir; 
                                        // for (int i = 0; i < g_bit.size(); i++){
                                        //         guide_correct_dir += convertBitsetToChar(g_bit[i]);
                                        //         target_correct_dir += convertBitsetToChar(t_bit[i]);
                                        //         if ((g_bit[i] & t_bit[i]) == 0 && ((!g_bit[i].none()) || (!t_bit[i].none()))){ //mismatch case and no '-' present
                                        //                 target_correct_dir[i] = target_correct_dir[i] + 32 ;
                                                                        
                                        //         }
                                        // }
                                      
                                        // vecInGuide.push_back(guide_correct_dir);
                                        // vecTargetOfGuide.push_back(target_correct_dir);
                                        vecInGuide.emplace_back(g);       // save guide
                                        vecTargetOfGuide.emplace_back(t); // save target
                                }

                        }else { //no target, only profile
                                if (bulType == 0)
                                        detailedOutputFast(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, profiling, ext_profiling, vecInGuide, vecTargetOfGuide);
                                else if (bulType > 0)
                                        detailedOutputFastBulgeDNA(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, bD, bulDNA, pd, pdm, vecInGuide, vecTargetOfGuide);
                                else
                                        detailedOutputFastBulgeRNA(guideI, g_bit, t_bit, g, t, bulType, mm, len_guide, bD, bulDNA, pr, prm, vecInGuide, vecTargetOfGuide);
                                
                        }
                        
                        index = (targetOnDNA[index].next + 1) * -1;


                        //Test for dup
                        // string target_correct_dir; 
                        // for (int i = 0; i < g_bit.size(); i++){
                        //         guide_correct_dir += convertBitsetToChar(g_bit[i]);
                        //         target_correct_dir += convertBitsetToChar(t_bit[i]);
                        //         if ((g_bit[i] & t_bit[i]) == 0 && ((!g_bit[i].none()) || (!t_bit[i].none()))){ //mismatch case and no '-' present
                        //                 target_correct_dir[i] = target_correct_dir[i] + 32 ;
                                                                        
                        //         }
                        // }

                        // duplicati << "\t" << target_correct_dir<<"\t" <<index << "\n";
                
                } while (index > -1); //duplicati << "\n";
        }
        else if (p[pos_tree].eqkid > 0)
                saveIndices( inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, d, bD, bR, bulType, targetOnDNA,  
                inGuide, targetOfGuide, vecTargetOfGuide, vecInGuide,bulgeType, directions,indices, mismatches, bulgeSize,
                guideI, profiling, ext_profiling, pd, pdm, pr, prm); //go to eqkid
}

// CONTROLLARE LA QUESTIONE BULGE CONTEMPORANEI IN DNA ED RNA
void nearsearch(char *s,int ti, int gi, vector<bitset<4>> &inGuide_bit, vector<bitset<4>> &targetOfGuide_bit, vector<Tnode> &p, int pos_tree, vector<bitset<4>> guide_bit, int pos_in_guide, int d, int bD, int bR, bool goToLoHi, int bulType,
                 Tleaf * t, char (&inGuide)[23], char (&targetOfGuide)[23], vector<string> &vTOG, vector<string> &vIG, vector<string> &bT, vector<char> &dirs, vector<int> &ind,  vector<int> &mismatches,  vector<int> &bulgeSize,
                 int guideI, vector<vector<vector<int>>> &profiling, vector<vector<vector<vector<vector<int>>>>> &ext_profiling, vector<vector<vector<int>>> &pd, vector<vector<vector<int>>> &pdm, vector<vector<vector<int>>>&pr, vector<vector<vector<int>>> &prm)
{               //char (&inGuide)[23], char (&targetOfGuide)[23]
        
        //cout << "guide: " << inGuide << ", target: " << targetOfGuide << endl;
        if (p[pos_tree].lokid > 0 && goToLoHi && (d > 0 || bD > 0 || bR > 0 || *s < p[pos_tree].splitchar)) // go to lokid // was guide_bit[pos_in_guide].to_ulong() < p[pos_tree].splitchar_bit.to_ulong()
                nearsearch(s, ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].lokid, guide_bit, pos_in_guide, d, bD, bR, true, bulType, t, 
                inGuide,targetOfGuide, vTOG, vIG, bT, dirs, ind,mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);

        if (p[pos_tree].hikid > 0 && goToLoHi && (d > 0 || bD > 0 || bR > 0 || *s > p[pos_tree].splitchar)) // go to hikid  //was guide_bit[pos_in_guide].to_ulong() > p[pos_tree].splitchar_bit.to_ulong()
                nearsearch(s, ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].hikid, guide_bit, pos_in_guide, d, bD, bR, true, bulType, t,
                inGuide,targetOfGuide, vTOG, vIG, bT, dirs, ind,mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);

        if (ti == (len_guide + bulDNA - bD))
        {            
                //string save1 = "save1";
               
                saveIndices( inGuide_bit, targetOfGuide_bit, p, pos_tree, d, bD, bR, bulType, t, inGuide,targetOfGuide, vTOG, vIG, bT, dirs, ind,mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm); 
        }
        else if (p[pos_tree].eqkid > 0)
        { // go to eqkid
                 
                inGuide[gi] = *s;                                        // save guide character
                inGuide_bit[gi] = guide_bit[pos_in_guide]; 
                gi++;
                if ((guide_bit[pos_in_guide] & p[pos_tree].splitchar_bit) != 0)
                { //MATCH CASE
                        targetOfGuide[ti] = p[pos_tree].splitchar;                     // save target character uppercase
                        targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit; 
                        ti++;
                        nearsearch(s+1,ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, guide_bit, pos_in_guide+1, d, bD, bR, true, bulType, t, 
                        inGuide,targetOfGuide, vTOG, vIG, bT, dirs, ind,mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);
                        //targetOfGuide[ti--] = 0;
                        //targetOfGuide_bit[ti--] = bitset<4>("0000");
                        ti--;
                }
                else if (d > 0)
                { // MISMATCH case
                        targetOfGuide[ti] = p[pos_tree].splitchar+32;                       // save target character lowercase
                        targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit; //remeber add -1
                        ti++;
                        nearsearch(s+1,ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, guide_bit, pos_in_guide+1, d - 1, bD, bR, true, bulType, t, 
                        inGuide,targetOfGuide, vTOG, vIG, bT, dirs, ind,mismatches, bulgeSize, guideI, profiling, ext_profiling,pd, pdm, pr, prm);
                        //targetOfGuide[ti--] = 0;
                        //targetOfGuide_bit[ti--] = bitset<4>("0000");
                        ti--;
                }
                //if (bR > 0 && *(s+1) && bulType < 1) {    // BULGE RNA case
                if (bR > 0 && (pos_in_guide + 1) < len_guide  && bulType < 1)
                { // BULGE RNA case
                        targetOfGuide[ti] = '-';                                       // update last target character with '-'
                        targetOfGuide_bit[ti] = bitset<4>("0000"); // il char '-'
                        ti++;
                        nearsearch(s + 1,ti, gi, inGuide_bit, targetOfGuide_bit, p, pos_tree, guide_bit, pos_in_guide+1, d, bD, bR - 1, false, bulType - 1, t, 
                        inGuide,targetOfGuide, vTOG, vIG, bT, dirs, ind,mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);
                        //targetOfGuide[ti--] = 0;
                        //targetOfGuide_bit[ti--] = bitset<4>("0000");
                        ti--;
                }
                if (bD > 0 && bulType > -1)
                { // BULGE DNA case
                        inGuide[gi-1] = '-';                                           // update last guide character with '-'
                        inGuide_bit[gi - 1] = bitset<4>("0000"); //il char '-'  //ho già inserito il char prima, lo aggiorno con '-'
                        targetOfGuide[ti] = p[pos_tree].splitchar;                             // save target character uppercase
                        targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit;
                        ti++;
                        nearsearch(s, ti, gi, inGuide_bit, targetOfGuide_bit, p, p[pos_tree].eqkid, guide_bit, pos_in_guide, d, bD - 1, bR, true, bulType + 1, t, 
                        inGuide,targetOfGuide, vTOG, vIG, bT, dirs, ind,mismatches, bulgeSize,guideI, profiling, ext_profiling, pd, pdm, pr, prm);
                        //targetOfGuide[ti--] = 0;
                        //targetOfGuide_bit[ti--] = bitset<4>("0000");
                        ti--;
                }
                //inGuide[gi--] = 0;
                
                //inGuide_bit[gi--] = bitset<4>("0000");
                gi--;
                                                                // current node is a leaf
        }
        
        //} else if (p->eqkid < 1 && *s) {
        
        //else if (p[pos_tree].eqkid < 1 && *s) //works
        else if (p[pos_tree].eqkid < 1 && pos_in_guide == (len_guide - 1))
        { //current node is a leaf
               
                inGuide[gi] = *s;                                                      // save guide character
                inGuide_bit[gi] = guide_bit[pos_in_guide]; //rememeber add -1
                gi++;
                targetOfGuide[ti] = *s == p[pos_tree].splitchar ? p[pos_tree].splitchar:p[pos_tree].splitchar+32;        // save target character       
               
                targetOfGuide_bit[ti] = p[pos_tree].splitchar_bit; //remebere add -1
                ti++;
                 //d = *s == p->splitchar ? d:d-1;        // update distance
                d = (guide_bit[pos_in_guide] & p[pos_tree].splitchar_bit) != 0 ? d : d - 1;
                if (d > -1 && ti == (len_guide + bulDNA - bD))
                { // save results
                        //string save2 = "save2";
                        
                        saveIndices( inGuide_bit, targetOfGuide_bit, p, pos_tree, d, bD, bR, bulType, t, inGuide,targetOfGuide, vTOG, vIG, bT, dirs, ind,mismatches, bulgeSize, guideI, profiling, ext_profiling, pd, pdm, pr, prm);
                }
                d = (guide_bit[pos_in_guide] & p[pos_tree].splitchar_bit) != 0 ? d : d + 1;
                //targetOfGuide[ti--] = 0;
                ti--;
                gi--;
                //inGuide[gi--] = 0;
                //targetOfGuide_bit[ti--] = bitset<4>("0000");
                //inGuide_bit[gi--] = bitset<4>("0000");
        }
}

void resetGlbVar(int numLeaves)
{
        double start, end; //for timing
        cout << "Reset time:\t";
        start = omp_get_wtime();

        //numNodes = 0;

        // for (int i = 0; i < numLeaves; i++)
        // {
        //         delete[] targetOnDNA[i].guideDNA_bit;
        //         delete[] targetOnDNA[i].guideDNA;
        // }

        // free(targetOnDNA); // array of target on DNA

        //numLeaves = 0;

        
        //pairNuc.clear(); pairNuc.resize(2);
        // pairNuc_bit.clear();
        // pairNuc_bit.resize(2);
        // flag = 1;
        //in = '0';
        kid = 0;
        //pNode = 0;
        //char inGuide[23]; // non servirà
        //inGuide_bit.clear(); inGuide_bit.resize(23);
        //char targetOfGuide[23]; // non servirà
        //targetOfGuide_bit.clear(); targetOfGuide_bit.resize(23);

        // vecTargetOfGuide.clear();
        // vecInGuide.clear();
        // bulgeType.clear();

        // vecTargetOfGuide_bit.clear(), vecInGuide_bit.clear(); //each target is a vector of bitset

        // directions.clear();
        // indices.clear();
        // mismatches.clear();
        // bulgeSize.clear();

        end = omp_get_wtime();
        cout << end - start << "\n";
}


void test_func(string &sss){
        sss[5] = 'g';
        cout << "sss: " << sss<< endl;
}

int main(int argc, char **argv)
{
        DIR *d;
        dirent *dir;
        vector<string> fileList;
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

        if (omp_get_max_threads() < num_thr)
                num_thr = omp_get_max_threads();

        getline(pamfile,line);
        transform(line.begin(), line.end(), line.begin(), ::toupper);// uppercase of the pam
        int delimiter = line.find(" ");
        string pam = line.substr(0, delimiter);
        int pamlimit = stoi(line.substr(delimiter, line.length() - 1)); //number that identifies the PAM length: NNNNNNNNNNNNNNNNNNNNNGG (3)
        
        int pamlen = pam.length(); //length of the total PAM: (NNNNNNNNNNNNNNNNNNNNNGG) is 23 
        pamRNA = pam.substr(pamlen - pamlimit, pamlimit);
        globalstart = omp_get_wtime(); // start global time

        cout << "Load Guides:\t";

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
                        guideRNA[numGuide][en] = '\0'; //was 20
                        
                }
                else
                {
                        guideRNA_s.push_back(line);
                        reverse(line.begin(), line.end());
                        guideRNA.push_back((char *)malloc(21 * sizeof(char)));
                        copy(line.begin(), line.end(), guideRNA[numGuide]); // save Guide
                        guideRNA[numGuide][en] = '\0'; // was 20
                }
                numGuide++;
        }

        //Transform loaded guides into bitset 
        for (int i = 0; i < guideRNA.size(); i++)
        {
                vector<bitset<4>> tmp(en);
                for (int j = 0; j < en; j++)
                {
                        tmp[j] = bitset<4>(convertCharToBitset(guideRNA[i][j]));
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
                        if (!(strcmp(".tmp", dir->d_name)))
                        {
                                continue;
                        }
                        fileList.push_back(dir->d_name);
                }
        }
        closedir(d);
        
        
        string gen_dir_s (genome_dir);
        gen_dir_s += "/.tmp";
        if (!(opendir(gen_dir_s.c_str())) )
                const int dir_err = mkdir(gen_dir_s.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
        //rmdir(gen_dir_s.c_str());
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
        if (argv[8] == resultwriting){
                fileResults.open(name_result + ".targets.txt", std::ios_base::binary); //out
                fileResults << "#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge Size\n";
                create_profile = false;
        }
        else if (argv[8] == profilewriting){
                fileprofiling.open(name_result + ".profile.xls", std::ios_base::out); 
                file_ext_profiling.open(name_result + ".extended_profile.xls", std::ios_base::out);
                file_profiling_dna.open(name_result + ".profile_dna.xls", std::ios_base::out);
                file_profiling_rna.open(name_result + ".profile_rna.xls", std::ios_base::out);

                create_target = false;
        } else{
                fileResults.open(name_result + ".targets.txt", std::ios_base::binary);//out
                fileprofiling.open(name_result + ".profile.xls", std::ios_base::out); 
                file_ext_profiling.open(name_result + ".extended_profile.xls", std::ios_base::out);
                file_profiling_dna.open(name_result + ".profile_dna.xls", std::ios_base::out);
                file_profiling_rna.open(name_result + ".profile_rna.xls", std::ios_base::out);
                fileResults << "#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge Size\n";
        }
        
        //Create matrix for profiling
        //vector<vector<int>> profiling(en  + mm + 1, vector<int>(numGuide, 0)); //was en + 2 + mm + 1 --> en = lunghezza guida, 2 = col ont offt, mm + 1 = numero di mismatch + colonna 0 mms
        //en = lunghezza guida (BP cells), mm = numero mismatch (1mm, 2mm ... cells), +1 = cella con 0 mm
        
        vector<vector<vector<int>>> profiling(en  + mm + 1, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));
        
        //Create matrix for extended profiling
        //vector<vector<vector<vector<int>>>> ext_profiling(numGuide, vector<vector<vector<int>>>(mm + 1, vector<vector<int>>(4, vector<int>(en, 0)))); 
        //1 dim = id of guide, 2 dim = select matrix with an x number of mms in the target, 3 dim = select nucleotide (acgt), 4 dim = select position
        
        vector<vector<vector<vector<vector<int>>>>> ext_profiling(numGuide, vector<vector<vector<vector<int>>>>(mm + 1, vector<vector<vector<int>>>(4, vector<vector<int>>(en, vector<int>(num_thr,0))))); 
        
        //Matrix for dna profiling
        // vector<vector<int>> profiling_dna_mm(en + bulDNA  + (mm + 1) * 2 , vector<int>(numGuide, 0));
        // vector<vector<int>> profiling_dna(en + bulDNA + (mm + 1) * 2, vector<int>(numGuide, 0));
        //the additional mm + 1 are for 

        vector<vector<vector<int>>> profiling_dna_mm(en + bulDNA  + (mm + 1) * 2 , vector<vector<int>>(numGuide, vector<int>(num_thr,0)));
        vector<vector<vector<int>>> profiling_dna(en + bulDNA  + (mm + 1) * 2 , vector<vector<int>>(numGuide, vector<int>(num_thr,0)));

        //Matrix for rna profiling
        // vector<vector<int>> profiling_rna_mm(en   + mm + 1 + mm + 1, vector<int>(numGuide, 0));
        // vector<vector<int>> profiling_rna(en  + mm + 1 + mm + 1, vector<int>(numGuide, 0));
        
        vector<vector<vector<int>>> profiling_rna_mm(en   + mm + 1 + mm + 1, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)    ));
        vector<vector<vector<int>>> profiling_rna(en   + mm + 1 + mm + 1, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)    ));

        len_guide = en;

        int numNodes;
        int numLeaves;
        Tleaf *targetOnDNA;
        

        vector<string> vecTargetOfGuide, vecInGuide, bulgeType;
        //vector<vector<bitset<4>>> vecTargetOfGuide_bit, vecInGuide_bit; //each target is a vector of bitset
        vector<char> directions;
        vector<int> indices, mismatches, bulgeSize;
        
        // string inGuide(en + pamRNA.size() + bulDNA,'x');
        // string targetOfGuide(en + pamRNA.size() + bulDNA,'x');
        char inGuide[23];
        char targetOfGuide[23];
        //string inGuide ;
        //inGuide =  "XXXXXXXXXXXXXX";
        //string targetOfGuide ;
        //targetOfGuide= "XXXXXXXXXXXXXXXXXXXXXX";
        //         cout << "inguide: " << typeid(inGuide).name()  << ", targetofGuide" << typeid(targetOfGuide).name()  << endl; 

        // test_func(inGuide);
        
        //INserire tutti i contatori
        int file;
        int jk;
        int i;
        //#PARALLEL SECTION
        //FOR ALL THE CHR
        #pragma omp parallel private(numNodes, numLeaves, targetOnDNA,  vecTargetOfGuide, vecInGuide, bulgeType, directions, indices, mismatches, bulgeSize, inGuide, targetOfGuide, file, jk, i) num_threads(num_thr)
        #pragma omp master
        for (file = 0; file < fileList.size(); file++)
        {
           #pragma omp task
           {
                string gen_dir(genome_dir);
                int last_under = fileList[file].find_last_of("_");
                int first_under = fileList[file].find_first_of("_");

                string chrName = fileList[file].substr(first_under + 1, last_under - first_under - 1);

                ifstream fileTree;
                // create tree
                vector<Tnode> albero;

                // read TST from file
                cout << "Load TST:\t";
                start = omp_get_wtime();
                
                targetOnDNA = loadTST(gen_dir + "/" + fileList[file], albero, fileTree, numNodes, numLeaves);
                //cout << targetOnDNA[0].guideIndex << endl;
                end = omp_get_wtime();
                cout << end - start << "\n";

                // Search guide and pam in the TST
                cout << "Nearsearch TST:\t";
                start = omp_get_wtime();

                double newguide = numGuide;
                int group_guide = ceil(newguide / 1000);

                for (jk = 0; jk < group_guide; jk++)
                {

                        int inizio = jk * 1000;
                        int fine = MIN(((jk + 1) * 1000), numGuide);
                        int mm2 = mm;
                        int bulDNA2 = bulDNA;
                        int bulRNA2 = bulRNA;
                        
                        //Tolto targetOnDna
                        vector<bitset<4>> inGuide_bit(en + pamRNA.size() + bulDNA); //mettere fuori e dichiararle private
                        vector<bitset<4>> targetOfGuide_bit(en + pamRNA.size() + bulDNA);
                        //#pragma omp for
                        for (i = inizio; i < fine; i++)
                        {
                                int ti = 0;
                                int gi = 0;
                                //#pragma omp single
                                {
                                        //#pragma omp task
                                        {
                                                nearsearch(guideRNA[i],ti, gi, inGuide_bit, targetOfGuide_bit, albero, 0, guideRNA_bit[i], 0, mm2, bulDNA2, bulRNA2, true, 0, 
                                                targetOnDNA, inGuide, targetOfGuide,vecTargetOfGuide, vecInGuide, bulgeType, directions, indices, mismatches, bulgeSize,
                                                i, profiling, ext_profiling, profiling_dna, profiling_dna_mm, profiling_rna, profiling_rna_mm); 
                                        }
                                }
                                //CONTROLLO
                                //inGuide_bit.clear(); inGuide_bit.resize(23);
                                //targetOfGuide_bit.clear(); targetOfGuide_bit.resize(23);
                        }
                        if (create_target && indices.size()>0){
                                string name_partial_res = gen_dir_s + "/"+ fileList[file] + "_" + to_string(omp_get_thread_num()) + "_" + to_string(jk) + ".txt";
                                ofstream partial_res(name_partial_res, std::ios_base::out);
                        
                                for (i = 0; i < indices.size(); i++)
                                {
                                        
                                        //string test = bulgeType[i] + "\t" + vecInGuide[i] +  "\t" + vecTargetOfGuide[i];
                                        //cout << "\n\n" << test << "\n\n";
                                        partial_res << bulgeType[i] << "\t" << vecInGuide[i] << "\t" << vecTargetOfGuide[i] << "\t" << chrName << "\t" << indices[i] << "\t" << directions[i] << "\t" << mismatches[i] << "\t" << bulgeSize[i] << "\n";
                                        
                                }
                        }
                       

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

                resetGlbVar(numLeaves);
                 for (i = 0; i < numLeaves; i++)
                {
                        delete[] targetOnDNA[i].guideDNA_bit;
                        delete[] targetOnDNA[i].guideDNA;
                }

                free(targetOnDNA);
           }
        }

        //Save profiling results
        if (create_profile){
                //Cols labels for profiling output file
                fileprofiling << "GUIDE\t";
                for (int i = 0; i < en; i++)
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

                //Cols labels for DNA profiling output file
                file_profiling_dna << "GUIDE\t";
                for (int i = 0; i < (en + bulDNA); i++)
                {
                        file_profiling_dna << "BP\t";
                }
                file_profiling_dna << "\t"
                              << "ONT\tOFFT\tOFFT/MM\t\t";
                for (int i = 0; i <= mm; i++)
                {
                        file_profiling_dna << i << "MM(1)\t";
                        file_profiling_dna << i << "MM(2)\t";
                }
                file_profiling_dna << "\n";
        
                //Cols labels for RNA profiling output file
                file_profiling_rna << "GUIDE\t";
                for (int i = 0; i < (en); i++)
                {
                        file_profiling_rna << "BP\t";
                }
                file_profiling_rna << "\t"
                              << "ONT\tOFFT\tOFFT/MM\t\t";
                for (int i = 0; i <= mm; i++)
                {
                        file_profiling_rna << i << "MM(1)\t";
                        file_profiling_rna << i << "MM(2)\t";
                }
                file_profiling_rna << "\n";
                for (int i = 0; i < numGuide; i++)
                {
                        saveProfileGuide(guideRNA_s[i], i, mm, en, bulDNA, profiling, ext_profiling, profiling_dna, profiling_dna_mm, profiling_rna, profiling_rna_mm,
                                        fileprofiling, file_ext_profiling, file_profiling_dna, file_profiling_rna, num_thr);
                }
                fileprofiling.close();
                file_ext_profiling.close();
                file_profiling_dna.close();
                file_profiling_rna.close();
        }
        
        if(create_target){
                //Get all file from tmp and joi to fileResults
                DIR *d1;
                dirent *dir1;
                vector<string> fileListTmp;
                d1 = opendir(gen_dir_s.c_str());
                if (d1)
                {
                        while ((dir1 = readdir(d1)) != NULL)
                        {
                                if (!(strcmp(".", dir1->d_name)))
                                {
                                        continue;
                                }
                                if (!(strcmp("..", dir1->d_name)))
                                {
                                        continue;
                                }
                                
                                fileListTmp.push_back(dir1->d_name);
                        }
                }
                closedir(d1);                
                for (int i = 0; i < fileListTmp.size(); i++){
                        
                        string res = gen_dir_s + "/" + fileListTmp[i]; 
                        //cout << "res " <<  res << endl;
                        ifstream a (res.c_str(), std::ios_base::binary);
                       
                                
                        fileResults << a.rdbuf() ;
                        
                        a.close();
                        remove (res.c_str());
                }

                
                fileResults.close();


        }
        //duplicati.close();

        //test join file
        // ofstream ia("a.txt",std::ios_base::out);
        // ia << "test\ttest\ttest\ntest_a\ttest_a\ntest_a2\ttest_2";
        // ia << "test\ttest\ttest\ntest_a\ttest_a\ntest_a2\ttest_2";
        // ofstream ib ("b.txt", std::ios_base::out);
        // ib << "ib ibb";
        // ofstream id("d.txt",std::ios_base::out);
        // id.close();
        // ib.close();
        // ia.close();
        // ifstream if_a("a.txt", std::ios_base::binary);
        // ifstream if_b("b.txt", std::ios_base::binary);
        // ifstream if_d("d.txt", std::ios_base::binary);
        // ofstream of_c("c.txt", std::ios_base::binary);
        
        // of_c << if_a.rdbuf() << if_d.rdbuf()<<if_b.rdbuf();
        
        // ifstream res ("test_output.targets.txt", std::ios_base::binary);
        // ofstream of_c("c.txt", std::ios_base::app);

        // for (int i = 0; i < 10; i++){
        //                 ifstream res ("test_output.targets.txt", std::ios_base::binary);

        //         of_c << res.rdbuf();
        // }
        
        cout << "\nC++ end" << endl;
        return 0;
}
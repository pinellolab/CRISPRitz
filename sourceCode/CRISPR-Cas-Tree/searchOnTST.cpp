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
#include "include/indexFormat.h"
#include <stdint.h>
#include <dirent.h> //for directory reading files
#include <array>
#include <sys/stat.h>
#include <iomanip>	  // std::setprecision
#include <sys/stat.h> //for file dimension

#include <chrono>
#include <unordered_set>  
#include <map>
#include <set>
#include <climits>
#include <iterator>

using namespace std;

/***************************************
 * -18/02/2020
 * 	-Patched PAM at beginning
 * 	-Saving results every 50 guides: ACTIVE
 * 	-Added column Cluster Position
 * 
 * Note: all databases created before 25/10/2019 are missing the length of the guide in the .bin file. This search is not compatible.
 * To avoid this issue, replace the read of indexGuideLen in loadTST
 * fileTree.read((char *)&indexGuideLen, sizeof(int));
 * with
 * indexGuideLen = 20;
 * (or the right size of the guide the tree was build), and skip the guide-length
 * check and note in main(), which read it from the same place.
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
	int guideIndex; // window start ws, decoded from either format by loadTST
	int next;
	unsigned char info;      // LEAF_* bits below, decoded by loadTST
	unsigned char pad_valid; // pad bases an alignment may consume (format 1: all)
	// (the two bytes sit in what was padding: a leaf stays 32 bytes, as in 2.9.1)
	char *guideDNA;
	bitset<4> *guideDNA_bit;
} Tleaf;

// Tleaf.info. Format 1 stores the strand in the sign of the leaf index, format 2 in
// a flag; loadTST turns both into these, so the emit path reads one form.
static const unsigned char LEAF_LEFT = 0x1;      // window = flank|guide|pad (PAM-left)
static const unsigned char LEAF_FLANK_OFF = 0x2; // part of the flank is off the contig
static const unsigned char LEAF_SKIP = 0x4;      // never emitted: a format-1 index-0 artefact (IDX-14), or a region hook's
static const unsigned char LEAF_REGION = 0x8;    // dropped by a region hook (a skip or keep BED, docs/INDEX_FORMAT.md section 13)

// What each thread needs to know about the .bin it is searching, set by loadTST:
// the index format, the max_bulges its coordinates are decoded with (format 2:
// its header's; format 1: the command line's, as 2.9.1 does), whether its build
// PAM is all N, and offset = the index's guide length minus the search's (2.9.1's
// offset_guide_len, a global then: each thread now keeps the one of its own .bin).
struct IndexFacts
{
	int format;
	int mb;
	bool pamless;
	int offset;
};
vector<IndexFacts> idxFacts;

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

vector<char *> guideRNA; // input guide RNA
// vector<string> guideRNA;				// input guide RNA
vector<string> guideRNA_s;				//input guide RNA as strings;
vector<vector<bitset<4>>> guideRNA_bit; //input guide RNA in bitset
int len_guide;
int max_bulges = 0;
int max_edits = -1; // total-edits cap (mismatches + DNA + RNA bulges), CRISPRme issue #107; -1 = disabled
bool pam_at_start;

// ---- On-the-fly PAM filter (opt-in) --------------------------------------
// When enabled, a PAM-less index (e.g. NNN) can serve a more restrictive
// search PAM (e.g. NGG/NRG/NRR/NGA/NAG/NNNN) WITHOUT a dedicated per-PAM index:
// at leaf emission, each candidate's stored concrete flanking PAM is compared
// against the IUPAC-expanded search PAM, and entries whose flank is incompatible
// with the search PAM are SKIPPED BEFORE they are written/scored. This is a
// leaf-time skip-before-write, NOT a tree-descent prune: the PAM-less tree is
// still traversed in full; the win is that excluded-PAM candidates are never
// emitted, scored or post-processed. The targets list is identical, as a set, to
// searching a dedicated index built for that restrictive PAM (see
// test/run_pam_filter_test.sh). The profile outputs are not: they count the
// engine's duplicate rows, whose number depends on the tree's shape, and their
// guide label carries the search PAM padded to the stored width K.
//
// Opt-in via: the 12th positional arg to searchTST (1 = on; crispritz.py
// `search --prune-pam`), OR the environment variable CRISPRITZ_PAM_PRUNE set to
// anything but empty or "0". Default (flag/env absent) is a strict no-op: the filter branch is never
// entered, so behaviour is byte-identical to upstream (see golden/minimal CI).
bool PAM_PRUNE = false;
int prune_pad = 0; // (K-L) distal N-pad width of a padded --prune-pam search;
                   // added to minus-strand coords (3' PAM) so output == a dedicated L-PAM index
vector<bitset<4>> searchPamBit; // IUPAC-expanded bitsets of the search PAM region (pamRNA)

// Local copy of the IUPAC->bitset conversion used BUILD-side in pamGenerator.cpp
// (same convention: A=0001 C=0010 G=0100 T=1000 N=1111 ...). Kept local so the
// searchTST link line is unchanged (pamGenerator.cpp is not linked into searchTST).
static vector<bitset<4>> pamPrune_bit_conversion(const string &PAM)
{
	vector<bitset<4>> pam_bit;
	for (size_t i = 0; i < PAM.length(); ++i)
	{
		char c = PAM[i];
		if (c == 'A')      pam_bit.push_back(bitset<4>(string("0001")));
		else if (c == 'C') pam_bit.push_back(bitset<4>(string("0010")));
		else if (c == 'G') pam_bit.push_back(bitset<4>(string("0100")));
		else if (c == 'T') pam_bit.push_back(bitset<4>(string("1000")));
		else if (c == 'N') pam_bit.push_back(bitset<4>(string("1111")));
		else if (c == 'R') pam_bit.push_back(bitset<4>(string("0101")));
		else if (c == 'Y') pam_bit.push_back(bitset<4>(string("1010")));
		else if (c == 'S') pam_bit.push_back(bitset<4>(string("0110")));
		else if (c == 'W') pam_bit.push_back(bitset<4>(string("1001")));
		else if (c == 'K') pam_bit.push_back(bitset<4>(string("1100")));
		else if (c == 'M') pam_bit.push_back(bitset<4>(string("0011")));
		else if (c == 'B') pam_bit.push_back(bitset<4>(string("1110")));
		else if (c == 'D') pam_bit.push_back(bitset<4>(string("1101")));
		else if (c == 'H') pam_bit.push_back(bitset<4>(string("1011")));
		else if (c == 'V') pam_bit.push_back(bitset<4>(string("0111")));
		else               pam_bit.push_back(bitset<4>(string("1111"))); // treat unknown as N
	}
	return pam_bit;
}

// Returns true if the stored leaf PAM (guideDNA_bit, indexed the SAME way as
// pamRNA/searchPamBit) is INCOMPATIBLE with the search PAM at any position, i.e.
// the entry should be pruned. Uses the BUILD-side convention: a position is a
// mismatch iff (stored & search).none().
static inline bool pamPrune_reject(const bitset<4> *leafPam)
{
	// Only the positions the search PAM actually covers are tested. The N-pad that
	// widens the search PAM to the stored flank width states no requirement, while a
	// stored 'N' is the 0000 nibble and intersects nothing -- so testing a padded
	// position would reject targets that a dedicated index of the real PAM width
	// keeps. Within the real width a stored 'N' still prunes, which is what a
	// dedicated index does when it matches its PAM on the genome.
	size_t first = pam_at_start ? (size_t)prune_pad : 0;
	size_t last = searchPamBit.size() - (pam_at_start ? 0 : (size_t)prune_pad);
	for (size_t j = first; j < last; ++j)
	{
		if ((leafPam[j] & searchPamBit[j]).none())
			return true; // incompatible at this PAM position -> prune
	}
	return false;
}
// ---- END PROTOTYPE globals ----

// The build PAM region, from the .bin file-name prefix (the only record of it in a
// format-1 index), as bitsets; and whether it constrains any position. Used to
// recognise the format-1 index-0 artefact (see loadTST).
vector<bitset<4>> buildPamBit;
bool buildPamConstrained = false;

// Whether the emit path must skip this leaf for an alignment that used dnaUsed DNA
// and rnaUsed RNA bulges:
//  - a format-1 index-0 artefact (LEAF_SKIP);
//  - --prune-pam: a flank incompatible with the search PAM (as before);
//  - no --prune-pam on a PAM-less index: a flank that runs off the contig. Such a
//    leaf exists only so a shorter PAM can be served (its real PAM bases are tested
//    above); without the filter it would be reported with a flank that is not there
//    and, on the minus strand, a coordinate before the contig start;
//  - an alignment that consumed more pad bases than lie on the contig before an N or
//    its end (pad_valid). It consumed len_guide + dnaUsed - rnaUsed stored bases, of
//    which the index's guide region holds the first len_guide + offset (offset = the
//    index's guide length minus the search's), so dnaUsed - rnaUsed - offset of them
//    are pad bases: a search guide longer than the index's reaches into the pad even
//    without a bulge, a shorter one only with DNA bulges. Format 1 stored only whole
//    windows, so it never skips there.
static inline bool leafSkipped(const Tleaf &leaf, int thr, int dnaUsed, int rnaUsed)
{
	if (leaf.info & LEAF_SKIP)
		return true;
	if (PAM_PRUNE)
	{
		if (pamPrune_reject(leaf.guideDNA_bit))
			return true;
	}
	else if (idxFacts[thr].pamless && (leaf.info & LEAF_FLANK_OFF))
		return true;
	return dnaUsed - rnaUsed - idxFacts[thr].offset > (int)leaf.pad_valid;
}

// ---- Region hooks: a skip BED, a keep BED, the leaf dump (docs/INDEX_FORMAT.md section 13) ----
// A search may be given a BED of intervals (contig, 0-based start, end exclusive) and
// SKIP every site whose core meets an interval (crispritz.py search --skip-bed, the 13th
// argument or CRISPRITZ_SKIP_BED), or KEEP only those sites (--keep-bed, the 14th argument
// or CRISPRITZ_KEEP_BED). The core of a site is the n = Gs + L bases of its bulge-free
// alignment: the L bases of the search PAM and the Gs search-guide bases next to them,
// whatever the index's own guide length G and flank width K. With off = G - Gs (section 1)
// and pp = K - L for a 3' PAM under --prune-pam (prune_pad; 0 otherwise), a leaf's core is
//   PAM-left  window [flank K][guide G][pad mb]:  [ws + pp,      ws + pp + n)
//   PAM-right window [pad mb][guide G][flank K]:  [ws + mb + off, ws + mb + off + n)
// which is the emit formulas of saveIndices run backwards: every line the leaf emits
// (Position p, d gaps in its crRNA, r in its DNA, n = its crRNA's length less d) has the
// core [p, p + n) if it is PAM-left and [p + d - r, p + d - r + n) if it is PAM-right. So
// a skip drops exactly the lines whose core, read from the line, meets an interval, and a
// keep keeps exactly those. A dropped leaf is marked LEAF_SKIP | LEAF_REGION when its tree
// is loaded, so the emit path drops it as it drops a format-1 index-0 artefact; and the
// traversal enters no subtree that holds no other leaf (two bits a node, computed once the
// tree is loaded; CRISPRITZ_SKIP_PRUNE=0 turns that off, which changes no output). The leaf
// dump (--leaf-dump, the 15th argument or CRISPRITZ_LEAF_DUMP) writes the stored letters of
// every leaf a skip BED drops. With no BED none of this runs: every check is behind
// regionActive, and the search is the one without the hooks, byte for byte.
typedef vector<pair<long long, long long>> RegionIntervals;
static map<string, RegionIntervals> regionIv; // per contig: sorted, merged, non-empty
static bool regionActive = false;             // a skip or a keep BED is in force
static bool regionKeep = false;               // the BED is a keep BED
static bool regionPruneOn = true;             // subtree pruning (CRISPRITZ_SKIP_PRUNE=0: off)
static bool regionPruneActive = false;        // regionActive && regionPruneOn
static string leafDumpPath;                   // empty: no leaf dump
static vector<char> regionTree;               // per thread: the loaded tree holds region-dropped leaves
static vector<vector<unsigned char>> nodeLive; // per thread, per node: bit 0 a live leaf below it (lo, hi
											   // or eq), bit 1 a live leaf through its eq kid
static vector<long long> nRegionSkipped, nRegionKept; // per thread
struct LeafDumpLine
{
	string contig;
	long long ws;
	char side; // 'L' PAM-left window, 'R' PAM-right
	string letters;
	bool operator<(const LeafDumpLine &o) const
	{
		if (contig != o.contig)
			return contig < o.contig;
		if (ws != o.ws)
			return ws < o.ws;
		if (side != o.side)
			return side < o.side;
		return letters < o.letters;
	}
};
static vector<vector<LeafDumpLine>> leafDump; // per thread

// The contig a tree file <PAM>_<contig>_<part>.bin holds (main reads it the same way).
static string binContig(const string &fileName)
{
	size_t first = fileName.find_first_of('_'), last = fileName.find_last_of('_');
	if (first == string::npos || last == string::npos || last <= first)
		return string();
	return fileName.substr(first + 1, last - first - 1);
}

// The argument at argv[pos] if given and not empty, else the environment variable env.
static string regionHookArg(int argc, char **argv, int pos, const char *env)
{
	if (argc > pos && argv[pos][0] != '\0')
		return string(argv[pos]);
	const char *v = getenv(env);
	return v != nullptr ? string(v) : string();
}

// Reads a region BED: tab-separated contig, start, end (0-based, end exclusive), more
// columns ignored, blank lines and lines starting with '#' skipped, a trailing CR
// dropped. A line that is not that, an end not after its start, or a contig that is
// not one of the index's (contigs) refuses the whole file: a BED read in part would
// skip, or keep, part of what it names, silently.
static bool loadRegionBed(const string &path, const set<string> &contigs, string &err, long long &nIntervals)
{
	ifstream f(path);
	if (!f)
	{
		err = "cannot be opened";
		return false;
	}
	string line;
	long long ln = 0;
	nIntervals = 0;
	while (getline(f, line))
	{
		ln++;
		if (!line.empty() && line[line.size() - 1] == '\r')
			line.erase(line.size() - 1);
		if (line.empty() || line[0] == '#')
			continue;
		vector<string> col;
		size_t a = 0;
		while (col.size() < 3)
		{
			size_t b = line.find('\t', a);
			col.push_back(line.substr(a, b == string::npos ? string::npos : b - a));
			if (b == string::npos)
				break;
			a = b + 1;
		}
		string where = "line " + to_string(ln) + ": ";
		if (col.size() < 3)
		{
			err = where + "fewer than 3 tab-separated columns (contig, start, end): '" + line.substr(0, 80) + "'";
			return false;
		}
		long long v[2];
		for (int k = 0; k < 2; k++)
		{
			const string &s = col[k + 1];
			if (s.empty() || s.size() > 18 || s.find_first_not_of("0123456789") != string::npos)
			{
				err = where + (k ? "the end" : "the start") + " '" + s.substr(0, 40) + "' is not a non-negative integer";
				return false;
			}
			v[k] = atoll(s.c_str());
		}
		if (col[0].empty())
		{
			err = where + "an empty contig name";
			return false;
		}
		if (v[1] <= v[0])
		{
			err = where + "the end " + to_string(v[1]) + " is not greater than the start " + to_string(v[0]);
			return false;
		}
		if (!contigs.count(col[0]))
		{
			string some;
			int k = 0;
			for (const string &c : contigs)
			{
				if (k++ == 5)
				{
					some += ", ...";
					break;
				}
				some += (k > 1 ? ", '" : "'") + c + "'";
			}
			err = where + "the contig '" + col[0].substr(0, 80) + "' is not one of the index's (" +
				  to_string(contigs.size()) + ": " + some + ")";
			return false;
		}
		regionIv[col[0]].push_back(make_pair(v[0], v[1]));
		nIntervals++;
	}
	if (f.bad())
	{
		err = "could not be read to its end";
		return false;
	}
	for (auto &kv : regionIv)
	{
		RegionIntervals &iv = kv.second, merged;
		sort(iv.begin(), iv.end());
		for (const auto &p : iv)
			if (!merged.empty() && p.first <= merged.back().second)
				merged.back().second = max(merged.back().second, p.second);
			else
				merged.push_back(p);
		iv.swap(merged);
	}
	return true;
}

// Whether [s, e) meets an interval of iv (sorted, merged).
static bool regionMeets(const RegionIntervals &iv, long long s, long long e)
{
	auto it = upper_bound(iv.begin(), iv.end(), make_pair(s, LLONG_MAX)); // the first interval starting after s
	if (it != iv.begin() && prev(it)->second > s)
		return true;
	return it != iv.end() && it->first < e;
}

// Marks the leaves of the tree just read (contig ctg) that the BED drops: in a skip, those
// whose core meets an interval; in a keep, the others. A format-1 index-0 artefact
// (already LEAF_SKIP, section 5) is no site and is never emitted: it is neither marked nor
// counted, so the leaves dropped are the sites dumped.
static void markRegionLeaves(int thr, const string &ctg, int numLeaves)
{
	auto f = regionIv.find(ctg);
	const RegionIntervals *iv = f != regionIv.end() ? &f->second : nullptr;
	const long long L = (long long)pamRNA.size() - prune_pad; // the search PAM's width
	const long long n = (long long)len_guide + L;
	const long long pp = pam_at_start ? 0 : prune_pad;
	const long long mbx = idxFacts[thr].mb, off = idxFacts[thr].offset;
	for (int i = 0; i < numLeaves; i++)
	{
		Tleaf &leaf = targetOnDNA[thr][i];
		if (leaf.info & LEAF_SKIP)
			continue;
		const long long ws = leaf.guideIndex;
		const long long cs = (leaf.info & LEAF_LEFT) ? ws + pp : ws + mbx + off;
		const bool meets = iv != nullptr && regionMeets(*iv, cs, cs + n);
		if (meets != regionKeep)
		{
			leaf.info |= LEAF_SKIP | LEAF_REGION;
			nRegionSkipped[thr]++;
			regionTree[thr] = 1;
		}
		else
			nRegionKept[thr]++;
	}
}

// Bit 0: a leaf that is not LEAF_SKIP lies below node n (through its lo, hi or eq kid);
// bit 1: one lies through its eq kid. Bottom-up, once a tree is loaded.
static unsigned char computeLive(vector<Tnode> &t, vector<unsigned char> &live, const vector<Tleaf> &leaves, int n)
{
	const Tnode &p = t[n];
	unsigned char eq = 0;
	if (p.eqkid < 0)
	{
		for (int idx = (p.eqkid + 1) * -1; idx >= 0; idx = (leaves[idx].next + 1) * -1)
			if (!(leaves[idx].info & LEAF_SKIP))
			{
				eq = 1;
				break;
			}
	}
	else if (p.eqkid > 0)
		eq = computeLive(t, live, leaves, p.eqkid) & 1;
	unsigned char any = eq;
	if (p.lokid > 0)
		any |= computeLive(t, live, leaves, p.lokid) & 1;
	if (p.hikid > 0)
		any |= computeLive(t, live, leaves, p.hikid) & 1;
	live[n] = any | (eq << 1);
	return live[n];
}

static char iupacComplement(char c)
{
	switch (c)
	{
	case 'A': return 'T';
	case 'T': return 'A';
	case 'C': return 'G';
	case 'G': return 'C';
	case 'R': return 'Y';
	case 'Y': return 'R';
	case 'K': return 'M';
	case 'M': return 'K';
	case 'B': return 'V';
	case 'V': return 'B';
	case 'D': return 'H';
	case 'H': return 'D';
	default: return c; // S, W, N, Q and '-' are their own complement
	}
}

// The leaf dump: walks the tree just loaded, rebuilding each string (the guide region
// and the pad, read from the PAM outward), and for every leaf a skip dropped writes its
// window on the + strand from ws, as the search would read it: the flank (in pamRNA
// order, the order loadTST decodes it in), the guide region and the pad. A pad base
// beyond pad_valid is written '-': the tree holds a placeholder there that no alignment
// may consume (section 3.2). So is a flank base off the contig: in a leaf flagged
// FLANK_OFF, the run of N (nibble 0) at the flank's distal end, which holds every flank
// base off the contig (section 3.3) and, at most, genomic N next to them, which match
// nothing. The other letters are the genome's, as the index stores them.
static void dumpLeaves(int thr, int node, int depth, string &path, const string &ctg)
{
	const Tnode &p = albero_glb[thr][node];
	if (p.lokid > 0)
		dumpLeaves(thr, p.lokid, depth, path, ctg);
	if (p.hikid > 0)
		dumpLeaves(thr, p.hikid, depth, path, ctg);
	if ((int)path.size() <= depth)
		path.resize(depth + 1);
	path[depth] = p.splitchar;
	if (p.eqkid > 0)
	{
		dumpLeaves(thr, p.eqkid, depth + 1, path, ctg);
		return;
	}
	for (int idx = (p.eqkid + 1) * -1; idx >= 0; idx = (targetOnDNA[thr][idx].next + 1) * -1)
	{
		const Tleaf &leaf = targetOnDNA[thr][idx];
		if (!(leaf.info & LEAF_REGION))
			continue;
		const int K = (int)pamRNA.size();
		string g(path, 0, depth + 1), fl(leaf.guideDNA, K), w;
		// in the guide's orientation: 5' PAM [flank][guide][pad], 3' PAM [pad][guide][flank]
		if (pam_at_start)
			w = fl + g;
		else
		{
			reverse(g.begin(), g.end());
			w = g + fl;
		}
		const int padInvalid = idxFacts[thr].mb - (int)leaf.pad_valid; // format 1: negative
		for (int z = 0; z < padInvalid && z < (int)g.size(); z++)
			w[pam_at_start ? w.size() - 1 - z : z] = '-';
		if (leaf.info & LEAF_FLANK_OFF)
			for (int z = 0; z < K; z++)
			{
				size_t at = pam_at_start ? (size_t)z : w.size() - 1 - z; // from the flank's distal end
				if (w[at] != 'N')
					break;
				w[at] = '-';
			}
		const bool minus = pam_at_start ? !(leaf.info & LEAF_LEFT) : (leaf.info & LEAF_LEFT);
		if (minus)
		{
			reverse(w.begin(), w.end());
			for (char &c : w)
				c = iupacComplement(c);
		}
		LeafDumpLine d;
		d.contig = ctg;
		d.ws = leaf.guideIndex;
		d.side = (leaf.info & LEAF_LEFT) ? 'L' : 'R';
		d.letters = w;
		leafDump[thr].push_back(d);
	}
}

// Once a tree is loaded: its liveness bits, and its leaf dump.
static void finishRegionTree(int thr, int numNodes, const string &ctg)
{
	if (!regionTree[thr])
		return;
	if (regionPruneOn)
	{
		nodeLive[thr].assign(numNodes, 0);
		computeLive(albero_glb[thr], nodeLive[thr], targetOnDNA[thr], 0);
	}
	if (!leafDumpPath.empty())
	{
		string path;
		dumpLeaves(thr, 0, 0, path, ctg);
	}
}
// ---- END region hooks ----

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

vector<vector<vector<int>>> profiling_both_mm;

vector<vector<string>> vecTargetOfGuide, vecInGuide, bulgeType;
vector<vector<char>> directions;
vector<vector<int>> indices, mismatches, bulgeSize, cluster_position;
vector<vector<string>> chrName_glb;




/**
 * This function reads all lines from the given file, keeps only the first occurrence
 * of each unique line (preserving the original order), and then rewrites the file 
 * without duplicates. The operation is performed in place: the input file is both
 * the source and the destination.
 *
 * @param path Path of the text file to clean from duplicate lines 
 *             (e.g. "results.targets.txt").
 */

void removeDuplicatesInPlace(const std::string &path){
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "Error opening file for reading: " << path << std::endl;
        return;
    }

    std::unordered_set<std::string> seen;
    std::vector<std::string> unique_lines;
    unique_lines.reserve(1024);

    std::string line;
    while (std::getline(in, line)) {
        if (seen.insert(line).second) { 
            unique_lines.push_back(std::move(line));
        }
    }
    in.close(); 

    std::ofstream out(path, std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "Error opening file for writing: " << path << std::endl;
        return;
    }
    for (const auto &l : unique_lines) out << l << '\n';
}

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
		unsigned char first = in & 0xF0; // keep only 4 left-most bits
		pairNuc_bit[0] = first >> 4;	 // shift them to the right to save them in a bitset<4>
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

// The tree of a file flagged CZI_FLAG_FOUR_ALLELE, read nibble by nibble. Its bytes are
// laid out as in any tree, but a slot's value depends on the slot (include/indexFormat.h):
// in the root, lo and hi slots 0 = no kid, 1..14 the IUPAC letters, 15 = Q; in an eq
// slot 15 = end of string, 1..14 the letters, 0 = Q. readPair cannot tell them apart: it
// reads every 1111 high nibble as an end of string and skips the byte's low nibble.
struct NibbleReader
{
	ifstream *f;
	unsigned char byte;
	bool have_low;
	int next()
	{
		if (have_low)
		{
			have_low = false;
			return byte & 0x0F;
		}
		char c;
		f->get(c);
		byte = (unsigned char)c;
		have_low = true;
		return byte >> 4;
	}
	void align() { have_low = false; } // an end of string in a high nibble ends its byte
};

// The IUPAC letter of a 4-bit A|C|G|T mask; 15 (all four) is Q, the fifth symbol.
static char nibbleSymbol(int n)
{
	static const char sym[16] = {'?', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'Q'};
	return sym[n & 0x0F];
}

// deSerialize for a file flagged CZI_FLAG_FOUR_ALLELE: n is the node's symbol, read from
// its parent's slot (an eq slot's 0 already turned into Q's 15).
static void deSerializeQ(NibbleReader &r, int n, int &pNode)
{
	int thr = omp_get_thread_num();
	int i = pNode;
	albero_glb[thr][i].splitchar = nibbleSymbol(n);
	albero_glb[thr][i].splitchar_bit = bitset<4>(n);
	int k = r.next(); // lo slot
	if (k != 0)
	{
		albero_glb[thr][i].lokid = ++pNode;
		deSerializeQ(r, k, pNode);
	}
	k = r.next(); // hi slot
	if (k != 0)
	{
		albero_glb[thr][i].hikid = ++pNode;
		deSerializeQ(r, k, pNode);
	}
	k = r.next(); // eq slot
	if (k == 0x0F)
	{
		r.align();
		r.f->read((char *)&(albero_glb[thr][i].eqkid), sizeof(int));
	}
	else
	{
		albero_glb[thr][i].eqkid = ++pNode;
		deSerializeQ(r, k == CZI_NIBBLE_Q_EQ ? CZI_NIBBLE_Q : k, pNode);
	}
}

// The header of a .bin (docs/INDEX_FORMAT.md). Format 1 (CRISPRitz <= 2.9.1) has
// none: its first int is the leaf count. Format 2 starts with CZI_MAGIC.
struct IndexHeader
{
	int format = 1;
	int version = 1;
	int headerLen = 0;
	int flags = 0;
	int guideLen = 0;
	int pamLen = 0;
	int mb = 0;
	string pam;
};

// Reads the header of the open .bin f, leaving f at the leaf count (the start of
// the file for format 1, header_len for format 2). Returns false, with the reason
// in err, for a file it must not read: a format it does not know, a malformed
// header, or a format-2 file that does not end in the end-of-index marker (cut
// short, e.g. by a full disk). The marker is checked before any of the tree is
// read, because deserializing a cut-short tree reads past its end.
static bool readIndexHeader(ifstream &f, IndexHeader &h, string &err)
{
	char m[4];
	if (!f.read(m, 4))
	{
		err = "is too short to be a CRISPRitz index (.bin) file";
		return false;
	}
	if (memcmp(m, CZI_MAGIC, 4) != 0)
	{
		h = IndexHeader();
		f.seekg(0);
		return true;
	}
	int v[6];
	if (!f.read((char *)v, sizeof(v)))
	{
		err = "has a cut-short index header";
		return false;
	}
	h.format = 2;
	h.version = v[0];
	h.headerLen = v[1];
	h.flags = v[2];
	h.guideLen = v[3];
	h.pamLen = v[4];
	h.mb = v[5];
	if (h.version > CZI_FORMAT_VERSION)
	{
		err = "is index format " + to_string(h.version) + ", which is newer than this CRISPRitz " +
			  "reads (formats 1 and 2): update CRISPRitz, or rebuild the index with this one";
		return false;
	}
	if (h.version != CZI_FORMAT_VERSION || (h.flags & ~CZI_FLAG_KNOWN) != 0 || h.pamLen < 1 ||
		h.pamLen > 1000 || h.headerLen != CZI_HEADER_FIXED + h.pamLen || h.guideLen < 1 ||
		h.mb < 0 || h.mb > CZI_MAX_BULGES)
	{
		err = "has a malformed index header (version " + to_string(h.version) + ", flags " +
			  to_string(h.flags) + ", PAM length " + to_string(h.pamLen) + ", header length " +
			  to_string(h.headerLen) + ")";
		return false;
	}
	h.pam.assign(h.pamLen, ' ');
	if (!f.read(&h.pam[0], h.pamLen))
	{
		err = "has a cut-short index header";
		return false;
	}
	char t[4];
	f.seekg(-4, ios::end);
	if (!f.read(t, 4) || memcmp(t, CZI_MAGIC, 4) != 0)
	{
		err = "is cut short: it does not end in the end-of-index marker";
		return false;
	}
	f.seekg(h.headerLen);
	return true;
}

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

	// main() has already read every header and refused a bad one; read this file's
	// again, for this thread. Format 2 decodes coordinates with its own max_bulges,
	// format 1 with the command line's, as 2.9.1 did.
	IndexHeader hdr;
	string hdrErr;
	if (!readIndexHeader(fileTree, hdr, hdrErr))
	{
		cerr << "ERROR: '" << path << "' " << hdrErr << "." << endl;
		exit(1);
	}
	idxFacts[thr].format = hdr.format;
	idxFacts[thr].mb = hdr.format == 2 ? hdr.mb : max_bulges;
	idxFacts[thr].pamless = hdr.format == 2 ? (hdr.flags & CZI_FLAG_PAMLESS) != 0 : !buildPamConstrained;
	// Q, the four-allele symbol, occurs only in a file flagged for it
	const bool quad = hdr.format == 2 && (hdr.flags & CZI_FLAG_FOUR_ALLELE) != 0;

	fileTree.read((char *)&numLeaves, sizeof(int)); // read number of leaves
	// Defense-in-depth: refuse to load a non-tree file (it should already be filtered to
	// *.bin at the directory scan). A garbage numLeaves would otherwise drive a huge
	// resize/alloc, and a later numNodes=0 would OOB-write in deSerialize.
	if (!fileTree || numLeaves <= 0)
	{
		cerr << "ERROR: '" << path << "' is not a valid CRISPRitz index (.bin) file "
			 << "(numLeaves=" << numLeaves << ")." << endl;
		exit(1);
	}
	targetOnDNA[thr].resize(numLeaves);

	// The guide length the index was built for. offset (2.9.1's offset_guide_len) is
	// added to the coordinates of PAM-right windows when the search guide is shorter
	// or longer than that, and counts the pad bases an alignment consumes (leafSkipped).
	// main() has checked it against a format-2 header.
	int indexGuideLen = 0;
	fileTree.read((char *)&indexGuideLen, sizeof(int));
	idxFacts[thr].offset = indexGuideLen - len_guide;

	char in = 0;
	for (int i = 0; i < numLeaves; i++)
	{																		 // fill array of targets on DNA
		fileTree.read((char *)&targetOnDNA[thr][i].guideIndex, sizeof(int)); // read index of target on DNA

		targetOnDNA[thr][i].guideDNA_bit = new bitset<4>[pamRNA.size()];
		targetOnDNA[thr][i].guideDNA = new char[pamRNA.size() + 1]; // initialize PAM size
		targetOnDNA[thr][i].guideDNA[pamRNA.size()] = '\0';
		unsigned char mask = 0;
		int k = 0;

		//read first char of PAM seq before entering the for cycle
		fileTree.get(in);

		for (int j = pamRNA.size() - 1; j > -1; j--)
		{
			if (k == 2) //when 2, one PAM char was read entirely
			{
				fileTree.get(in);
				k = 0;
			}
			
			//if PAM size is uneven, necessary to shift the last char to read the last nt of PAM seq
			if (j==0 && (pamRNA.size()%2))
				in <<= 4;

			//mask to read each time the correct 4bits of the char containing the PAM nt
			mask = in & 0xF0;//mask with 11110000
			in <<= 4;//shift to read the correct 4 bits

			switch (mask)
			{
			case 0x0:
				targetOnDNA[thr][i].guideDNA[j] = 'N';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("0000"); //intersects nothing
				break;
			case 0x10:
				targetOnDNA[thr][i].guideDNA[j] = 'A';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("0001");
				// cout << "A" <<endl;
				break;
			case 0x20:
				targetOnDNA[thr][i].guideDNA[j] = 'C';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("0010");
				// cout << "C"<<endl;
				break;
			case 0x40:
				targetOnDNA[thr][i].guideDNA[j] = 'G';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("0100");
				// cout << "G" <<endl;
				break;
			case 0x80:
				targetOnDNA[thr][i].guideDNA[j] = 'T';
				targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("1000");
				// cout << "T" <<endl;
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
			case 0xF0: //Q, all four bases: only in a file flagged CZI_FLAG_FOUR_ALLELE
				if (quad)
				{
					targetOnDNA[thr][i].guideDNA[j] = 'Q';
					targetOnDNA[thr][i].guideDNA_bit[j] = bitset<4>("1111");
					break;
				}
				cerr << "Error on reading the PAM, input char is not an IUPAC nucleotide: " << mask << endl;
				break;
			default:
				cerr << "Error on reading the PAM, input char is not an IUPAC nucleotide: " << mask << endl;
				break;
			}
			k++;
		}


		// std::string str(targetOnDNA[thr][i].guideDNA);
		// cout << "pam rna read " << str << endl;

		Tleaf &leaf = targetOnDNA[thr][i];
		fileTree.get(in); // format 1: '0' or '_' and the next leaf; format 2: the leaf flags
		if (hdr.format == 1)
		{
			if (in == '0')
				leaf.next = 0;
			else
				fileTree.read((char *)&leaf.next, sizeof(int));

			// The sign of a format-1 leaf index is the strand: < 0 is a PAM-left window
			// (3' PAM minus strand, 5' PAM plus strand) at -index, > 0 a PAM-right one.
			// Index 0 is ambiguous: the 2.9.1 builder makes it PAM-left for a 3' PAM and
			// PAM-right for a 5' one, so read it the same way (IDX-14; 2.9.1 read index
			// 0 as PAM-right for both, so a reverse-strand PAM at a contig's first base
			// came out on the + strand).
			int gi = leaf.guideIndex;
			leaf.guideIndex = gi < 0 ? -gi : gi;
			leaf.info = (gi < 0 || (gi == 0 && !pam_at_start)) ? LEAF_LEFT : 0;
			leaf.pad_valid = 255; // format 1 kept only windows that fit: the whole pad is real
			// That builder also took a window of the OTHER shape that starts at 0 (a 3'
			// PAM-right window whose guide starts at max_bulges, a 5' PAM-left one whose
			// PAM is at base 0) for index 0, and stored it from the wrong bases. Its flank
			// then matches the build PAM only where those bases hold a real PAM of the
			// shape index 0 is read as -- and then it is that real window, which the scan
			// found too. So an index-0 flank that fails the build PAM is such an artefact:
			// drop it. (An all-N build PAM: every window is real, nothing to drop.)
			if (gi == 0 && buildPamConstrained)
				for (size_t j = 0; j < buildPamBit.size() && j < pamRNA.size(); j++)
					if ((leaf.guideDNA_bit[j] & buildPamBit[j]).none())
					{
						leaf.info |= LEAF_SKIP;
						break;
					}
		}
		else
		{
			unsigned char lflags = (unsigned char)in;
			if (lflags & CZI_LEAF_HAS_NEXT)
				fileTree.read((char *)&leaf.next, sizeof(int));
			else
				leaf.next = 0;
			// format 2 stores the window start itself, signed, and the shape in a flag
			leaf.info = ((lflags & CZI_LEAF_PAM_LEFT) ? LEAF_LEFT : 0) |
						((lflags & CZI_LEAF_FLANK_OFF) ? LEAF_FLANK_OFF : 0);
			leaf.pad_valid = lflags & CZI_LEAF_PAD_MASK;
		}
	}

	// region hooks (a skip or keep BED): mark the leaves the BED drops
	string regionCtg;
	if (regionActive)
	{
		regionCtg = binContig(path.substr(path.find_last_of('/') + 1));
		regionTree[thr] = 0;
		markRegionLeaves(thr, regionCtg, numLeaves);
	}

	fileTree.read((char *)&numNodes, sizeof(int)); // read number of nodes
	if (!fileTree || numNodes <= 0)
	{
		cerr << "ERROR: '" << path << "' has an invalid TST node count (numNodes="
			 << numNodes << "); refusing to deserialize." << endl;
		exit(1);
	}
	vector<bitset<4>> pairNuc_bit(2);
	char pairNuc[2];
	uint8_t flag = 1;

	albero_glb[omp_get_thread_num()].resize(numNodes);
	int pNode = 0;
	if (quad)
	{
		NibbleReader r = {&fileTree, 0, false};
		deSerializeQ(r, r.next(), pNode); // the root's slot is no eq slot
		fileTree.close();
		if (regionActive)
			finishRegionTree(thr, numNodes, regionCtg);
		return;
	}
	fileTree.get(in);

	readPair(pairNuc_bit, pairNuc, in);

	flag = 0;
	deSerialize(albero, fileTree, pairNuc_bit, pairNuc, in, flag, pNode); // deserialize TST

	fileTree.close();
	if (regionActive)
		finishRegionTree(thr, numNodes, regionCtg);
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
 * @param bulTypeRNA Integer value that represents the type of the bulge: <0 is a RNA bulge, >0 is a DNA bulge, =0 means there are no bulges
 * (" X ") 
 * @param inGuide Array of char that contains the Guide in input and possible bulges
 * @param targetOfGuide Array of char that contains the Target found and possible bulges
 */
void saveIndices(Tnode *p, int d, int bD, int bR, int bulTypeRNA, int bulTypeDNA)
{
	if (p->lokid > 0)
		saveIndices(&albero_glb[omp_get_thread_num()][p->lokid], d, bD, bR, bulTypeRNA, bulTypeDNA); // go to lokid
	if (p->hikid > 0)
		saveIndices(&albero_glb[omp_get_thread_num()][p->hikid], d, bD, bR, bulTypeRNA, bulTypeDNA); // go to hikid

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
			for (int i = 0; i < (int)pamRNA.size() - prune_pad; i++)
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
			// On-the-fly PAM filter (mismatch path): skip this leaf entry BEFORE
			// emitting if its stored concrete flank PAM is incompatible with the
			// search PAM. Leaf-time skip-before-write (not a during-descent prune):
			// the index is PAM-less so the descent still visits every node; we only
			// avoid writing the excluded-PAM entries. Advancing `index` to the next
			// entry then `continue` mirrors the normal end-of-loop advance below.
			// leafSkipped also drops the leaves a format-2 index stores for boundary
			// sites when this alignment needs bases they lack (pad or flank off the
			// contig), and the format-1 index-0 artefact.
			if (leafSkipped(targetOnDNA[thr][index], thr, bulDNA - bD, bulRNA - bR))
			{
				index = (targetOnDNA[thr][index].next + 1) * -1;
				continue;
			}

			string t(targetOfGuide[thr]);     // take all the string
			// cout << "t subito assignment " << t << endl;

			t = t.substr(0, len_guide + bulDNA - bD); // take a substring of the string 
			// cout << "t substring " << t << endl;

			if (pam_at_start)
			{
				string tmp(targetOnDNA[thr][index].guideDNA);
				t.insert(0, tmp);
			}
			else
			{
				reverse(t.begin(), t.end());
				t.append(targetOnDNA[thr][index].guideDNA, (size_t)((int)pamRNA.size() - prune_pad));
			}

			// cout << "t pam assignment " << t << endl;

			vector<bitset<4>> t_bit(len_guide + bulDNA);
			int i = 0;
			for (i = 0; i < (len_guide + bulDNA - bD); i++)
			{
				t_bit[i] = targetOfGuide_bit[thr][i];
			}

			if (!pam_at_start)
			{
				reverse(t_bit.begin(), t_bit.begin() + len_guide + bulDNA - bD);
			}


			//Check what type of result to create
			if (create_target)
			{

				mismatches[thr].emplace_back(mm - d); // save mismatches
				if (bulTypeRNA == 0 && bulTypeDNA == 0)
				{ // NO BULGE case
					bulgeType[thr].emplace_back("X");
					bulgeSize[thr].emplace_back(0);
				}
				else
				{ // BULGE case
					if (bulTypeRNA != 0 && bulTypeDNA != 0){
						bulgeType[thr].emplace_back("RNA,DNA");
					}
					else
					{
						bulgeType[thr].emplace_back(bulTypeRNA > 0 ? "RNA" : "DNA");
					}

					bulgeSize[thr].emplace_back((bulRNA - bR) + (bulDNA - bD));
				}
				if (targetOnDNA[thr][index].info & LEAF_LEFT)
				{ // PAM-left window (3' PAM: minus strand; 5' PAM: plus strand)
					indices[thr].emplace_back(targetOnDNA[thr][index].guideIndex + (pam_at_start ? 0 : prune_pad));
					cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + (pam_at_start ? 0 : prune_pad));
					if (pam_at_start)
					{
						directions[thr].emplace_back('+');
					}
					else
					{
						directions[thr].emplace_back('-');
					}
				}
				else
				{ // PAM-right window (3' PAM: plus strand; 5' PAM: minus strand)
					indices[thr].emplace_back(targetOnDNA[thr][index].guideIndex + idxFacts[thr].mb - (bulDNA - bD) + (bulRNA - bR) + idxFacts[thr].offset);
					if (bulTypeRNA == 0 && bulTypeDNA == 0)
					{ // X case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + idxFacts[thr].mb - (bulDNA - bD) + (bulRNA - bR) + idxFacts[thr].offset);
					}
					else if (bulTypeRNA > 0 && bulTypeDNA == 0)
					{ //RNA case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + idxFacts[thr].mb - (bulDNA - bD) + (bulRNA - bR) + idxFacts[thr].offset - (bulRNA - bR));
					}
					else if (bulTypeRNA == 0 && bulTypeDNA > 0)
					{ //DNA case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + idxFacts[thr].mb - (bulDNA - bD) + (bulRNA - bR) + idxFacts[thr].offset + (bulDNA - bD));
					}
					else
					{ //DNA case and RNA case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + idxFacts[thr].mb - (bulDNA - bD) + (bulRNA - bR) + idxFacts[thr].offset + (bulDNA - bD) - (bulRNA - bR));
					}

					if (pam_at_start)
					{
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
					if (bulTypeRNA == 0 && bulTypeDNA == 0){
						//printf("\nno bulge\n");
						detailedOutputFast(guideI[thr], g_bit, t_bit, g, t, mm, len_guide, profiling, ext_profiling, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start);}
					else if (bulTypeRNA == 0 && bulTypeDNA > 0){
						//printf("\nbulge DNA\n");
						detailedOutputFastBulgeDNA(guideI[thr], g_bit, t_bit, g, t, mm, max_bulges, len_guide, bD, bulDNA, profiling_dna, profiling_dna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling);}
					else if (bulTypeRNA > 0 && bulTypeDNA == 0){
						//printf("\nbulge RNA\n");
						detailedOutputFastBulgeRNA(guideI[thr], g_bit, t_bit, g, t, mm, max_bulges, len_guide, bD, bulDNA, profiling_rna, profiling_rna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_rna, ext_profiling);}
					else{
						//printf("\nboth RNA and DNA bulges\n");
						detailedOutputFastBulgeBoth( guideI[thr], g_bit, t_bit, g, t, mm, max_bulges,len_guide,bD,bulDNA,profiling_dna,profiling_dna_mm,profiling_rna,profiling_rna_mm,
							profiling_both_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling_rna, ext_profiling);
						}
				}
				else
				{
					vecInGuide[thr].emplace_back(g);	   // save guide
					vecTargetOfGuide[thr].emplace_back(t); // save target
														   // cout << "target found " << t << endl;
				}
			}
			else
			{ //no target, only profile  !!!TO EDIT!!!!
			
				if (bulTypeRNA == 0 && bulTypeDNA == 0){
					//printf("\nno bulge\n");
					detailedOutputFast(guideI[thr], g_bit, t_bit, g, t, mm, len_guide, profiling, ext_profiling, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start);}
				else if (bulTypeRNA == 0 && bulTypeDNA > 0){
					//printf("\nbulge DNA\n");
					detailedOutputFastBulgeDNA(guideI[thr], g_bit, t_bit, g, t, mm, max_bulges, len_guide, bD, bulDNA, profiling_dna, profiling_dna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling);}
				else if (bulTypeRNA > 0 && bulTypeDNA == 0){
					printf("\nbulge RNA\n");
					detailedOutputFastBulgeRNA(guideI[thr], g_bit, t_bit, g, t, mm, max_bulges, len_guide, bD, bulDNA, profiling_rna, profiling_rna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_rna, ext_profiling);}
				else{
					//printf("\nboth RNA and DNA bulges\n");
					detailedOutputFastBulgeBoth( guideI[thr], g_bit, t_bit, g, t, mm, max_bulges,len_guide,bD,bulDNA,profiling_dna,profiling_dna_mm,profiling_rna,profiling_rna_mm,
						profiling_both_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling_rna, ext_profiling);
					}
			}

			index = (targetOnDNA[thr][index].next + 1) * -1;

		} while (index > -1);
	}
	else if (p->eqkid > 0)
		saveIndices(&albero_glb[omp_get_thread_num()][p->eqkid], d, bD, bR, bulTypeRNA, bulTypeDNA); //go to eqkid
}

void saveIndices2(Tnode *p, int d, int bD, int bR, int bulTypeRNA, int bulTypeDNA)
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
			for (int i = 0; i < (int)pamRNA.size() - prune_pad; i++)
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
			// On-the-fly PAM filter (bulge path): skip this leaf entry BEFORE
			// emitting if its stored concrete flank PAM is incompatible with the
			// search PAM. Same semantics as saveIndices; advance to the next entry.
			// leafSkipped also drops the leaves a format-2 index stores for boundary
			// sites when this alignment needs bases they lack (pad or flank off the
			// contig), and the format-1 index-0 artefact.
			if (leafSkipped(targetOnDNA[thr][index], thr, bulDNA - bD, bulRNA - bR))
			{
				index = (targetOnDNA[thr][index].next + 1) * -1;
				continue;
			}

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
				t.append(targetOnDNA[thr][index].guideDNA, (size_t)((int)pamRNA.size() - prune_pad));
			}

			vector<bitset<4>> t_bit(len_guide + bulDNA);
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
				if (bulTypeRNA == 0 && bulTypeDNA == 0)
				{ // NO BULGE case
					bulgeType[thr].emplace_back("X");
					bulgeSize[thr].emplace_back(0);
				}
				else
				{ // BULGE case
					if (bulTypeRNA != 0 && bulTypeDNA != 0){
						bulgeType[thr].emplace_back("RNA,DNA");
					}
					else
					{
						bulgeType[thr].emplace_back(bulTypeRNA > 0 ? "RNA" : "DNA");
					}

					bulgeSize[thr].emplace_back((bulRNA - bR) + (bulDNA - bD));
				}

				if (targetOnDNA[thr][index].info & LEAF_LEFT)
				{ // PAM-left window (3' PAM: minus strand; 5' PAM: plus strand)
					indices[thr].emplace_back(targetOnDNA[thr][index].guideIndex + (pam_at_start ? 0 : prune_pad));
					cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + (pam_at_start ? 0 : prune_pad));
					if (pam_at_start)
					{
						directions[thr].emplace_back('+');
					}
					else
					{
						directions[thr].emplace_back('-');
					}
				}
				else
				{ // PAM-right window (3' PAM: plus strand; 5' PAM: minus strand)
					indices[thr].emplace_back(targetOnDNA[thr][index].guideIndex + idxFacts[thr].mb - (bulDNA - bD) + (bulRNA - bR) + idxFacts[thr].offset);
					if (bulTypeRNA == 0 && bulTypeDNA == 0)
					{ // X case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + idxFacts[thr].mb - (bulDNA - bD) + (bulRNA - bR) + idxFacts[thr].offset);
					}
					else if (bulTypeRNA > 0 && bulTypeDNA == 0)
					{ //RNA case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + idxFacts[thr].mb - (bulDNA - bD) + (bulRNA - bR) + idxFacts[thr].offset - (bulRNA - bR));
					}
					else if (bulTypeRNA == 0 && bulTypeDNA > 0)
					{ //DNA case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + idxFacts[thr].mb - (bulDNA - bD) + (bulRNA - bR) + idxFacts[thr].offset + (bulDNA - bD));
					}
					else
					{ //DNA case and RNA case
						cluster_position[thr].emplace_back(targetOnDNA[thr][index].guideIndex + idxFacts[thr].mb - (bulDNA - bD) + (bulRNA - bR) + idxFacts[thr].offset + (bulDNA - bD) - (bulRNA - bR));
					}

					if (pam_at_start)
					{
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
				if (bulTypeRNA == 0 && bulTypeDNA == 0){
					//printf("\nno bulge\n");
					detailedOutputFast(guideI[thr], g_bit, t_bit, g, t, mm, len_guide, profiling, ext_profiling, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start);}
				else if (bulTypeRNA == 0 && bulTypeDNA > 0){
					//printf("\nbulge DNA\n");
					detailedOutputFastBulgeDNA(guideI[thr], g_bit, t_bit, g, t, mm, max_bulges, len_guide, bD, bulDNA, profiling_dna, profiling_dna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling);}
				else if (bulTypeRNA > 0 && bulTypeDNA == 0){
					//printf("\nbulge RNA\n");
					detailedOutputFastBulgeRNA(guideI[thr], g_bit, t_bit, g, t, mm, max_bulges, len_guide, bD, bulDNA, profiling_rna, profiling_rna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_rna, ext_profiling);}
				else{
					//printf("\nboth RNA and DNA bulges\n");
					detailedOutputFastBulgeBoth( guideI[thr], g_bit, t_bit, g, t, mm, max_bulges,len_guide,bD,bulDNA,profiling_dna,profiling_dna_mm,profiling_rna,profiling_rna_mm,
						profiling_both_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling_rna, ext_profiling);
					}
				}
				else
				{
					vecInGuide[thr].emplace_back(g);	   // save guide
					vecTargetOfGuide[thr].emplace_back(t); // save target
														   // cout << "target found " << t << endl;
				}
			}
			else
			{ //no target, only profile
				if (bulTypeRNA == 0 && bulTypeDNA == 0){
					//printf("\nno bulge\n");
					detailedOutputFast(guideI[thr], g_bit, t_bit, g, t, mm, len_guide, profiling, ext_profiling, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start);}
				else if (bulTypeRNA == 0 && bulTypeDNA > 0){
					//printf("\nbulge DNA\n");
					detailedOutputFastBulgeDNA(guideI[thr], g_bit, t_bit, g, t, mm, max_bulges, len_guide, bD, bulDNA, profiling_dna, profiling_dna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling);}
				else if (bulTypeRNA > 0 && bulTypeDNA == 0){
					//printf("\nbulge RNA\n");
					detailedOutputFastBulgeRNA(guideI[thr], g_bit, t_bit, g, t, mm, max_bulges, len_guide, bD, bulDNA, profiling_rna, profiling_rna_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_rna, ext_profiling);}
				else{
					//printf("\nboth RNA and DNA bulges\n");
					detailedOutputFastBulgeBoth( guideI[thr], g_bit, t_bit, g, t, mm, max_bulges,len_guide,bD,bulDNA,profiling_dna,profiling_dna_mm,profiling_rna,profiling_rna_mm,
						profiling_both_mm, vecInGuide[thr], vecTargetOfGuide[thr], pam_at_start, ext_profiling_dna, ext_profiling_rna, ext_profiling);
				}
			}

			index = (targetOnDNA[thr][index].next + 1) * -1;

		} while (index > -1);
	}
	else if (p->eqkid > 0)
		saveIndices2(&albero_glb[omp_get_thread_num()][p->eqkid], d, bD, bR, bulTypeRNA, bulTypeDNA); //go to eqkid
}

/**
 * Function that, given a Guide in input, searches it inside the TST. At first it goes recursively to the left child, then to the right child
 * and then compares the node character with the current Guide character.
 * 
 */
void nearsearch(char *s, Tnode *p, int pos_in_guide, int d, int bD, int bR, bool goToLoHi, int bulTypeRNA, int bulTypeDNA,
				const int thr)
{
	// CRISPRme issue #107: prune whole subtrees whose path already spends more than
	// the total-edits budget (mismatches + bulges). Edits only grow going down, and
	// lo/hi siblings share the same (d,bD,bR), so pruning here is safe and complete.
	if (max_edits >= 0 && (mm - d) + (bulDNA - bD) + (bulRNA - bR) > max_edits)
		return;
	// region hooks: enter no subtree without a leaf the BED keeps (bit 0), and no eq
	// branch without one (bit 1); the leaves themselves are dropped at emission
	unsigned char live = 3;
	if (regionPruneActive && regionTree[thr])
	{
		live = nodeLive[thr][p - albero_glb[thr].data()];
		if (!(live & 1))
			return;
	}
	if (p->lokid > 0 && goToLoHi)
	{ // go to lokid
		nearsearch(s, &albero_glb[thr][p->lokid], pos_in_guide, d, bD, bR, true, bulTypeRNA, bulTypeDNA, thr);
	}

	if (p->hikid > 0 && goToLoHi)
	{ // go to hikid
		nearsearch(s, &albero_glb[thr][p->hikid], pos_in_guide, d, bD, bR, true, bulTypeRNA, bulTypeDNA, thr);
	}

	if (ti[thr] == (len_guide + bulDNA - bD) &&
		(max_edits < 0 || (mm - d) + (bulDNA - bD) + (bulRNA - bR) <= max_edits))
	{

		save_function_start = omp_get_wtime();
		saveIndices(p, d, bD, bR, bulTypeRNA, bulTypeDNA);
		save_function_end = omp_get_wtime();
		save_function += save_function_end - save_function_start;
	}
	else if (p->eqkid > 0 && (live & 2))
	{ // go to eqkid

		inGuide[thr][gi[thr]] = *s; // save guide character
		inGuide_bit[thr][gi[thr]] = guideRNA_bit[guideI[thr]][pos_in_guide];  // guideRNA_bit = guideRNA in bitset
		gi[thr]++;
		if ((guideRNA_bit[guideI[thr]][pos_in_guide] & p->splitchar_bit) != 0)
		{												//MATCH CASE
			targetOfGuide[thr][ti[thr]] = p->splitchar; // save target character uppercase
			targetOfGuide_bit[thr][ti[thr]] = p->splitchar_bit;
			ti[thr]++;
			nearsearch(s + 1, &albero_glb[thr][p->eqkid], pos_in_guide + 1, d, bD, bR, true, bulTypeRNA, bulTypeDNA, thr);
			ti[thr]--;
		}
		else if (d > 0)
		{													 // MISMATCH case
			targetOfGuide[thr][ti[thr]] = p->splitchar + 32; // save target character lowercase
			targetOfGuide_bit[thr][ti[thr]] = p->splitchar_bit;
			ti[thr]++;
			nearsearch(s + 1, &albero_glb[thr][p->eqkid], pos_in_guide + 1, d - 1, bD, bR, true, bulTypeRNA, bulTypeDNA, thr);
			ti[thr]--;
		}

		if (bR > 0)            // bulge in the guide 
		{														 // BULGE RNA case
			targetOfGuide[thr][ti[thr]] = '-';					 // update last target character with '-'
			targetOfGuide_bit[thr][ti[thr]] = bitset<4>("0000"); // char '-'
			ti[thr]++;
			nearsearch(s + 1, p, pos_in_guide + 1, d, bD, bR - 1, true, bulTypeRNA+1 , bulTypeDNA, thr);
			ti[thr]--;
		}
		if (bD > 0)
		{													   // BULGE DNA case
			inGuide[thr][gi[thr] - 1] = '-';				   // update last guide character with '-'
			inGuide_bit[thr][gi[thr] - 1] = bitset<4>("0000"); //char '-'  //i've alredy inserted a char into inGuide, so it's updated with '-'
			targetOfGuide[thr][ti[thr]] = p->splitchar;		   // save target character uppercase
			targetOfGuide_bit[thr][ti[thr]] = p->splitchar_bit;
			ti[thr]++;
			nearsearch(s, &albero_glb[thr][p->eqkid], pos_in_guide, d, bD - 1, bR, true, bulTypeRNA, bulTypeDNA+1, thr);

			ti[thr]--;
		}
		gi[thr]--;
	}
	else if (p->eqkid < 1 && pos_in_guide == (len_guide - 1)) // 
	{ //current node is a leaf

		inGuide[thr][gi[thr]] = *s; // save guide character
		inGuide_bit[thr][gi[thr]] = guideRNA_bit[guideI[thr]][pos_in_guide];
		gi[thr]++;
		targetOfGuide[thr][ti[thr]] = (guideRNA_bit[guideI[thr]][pos_in_guide] & p->splitchar_bit) != 0 ? p->splitchar : p->splitchar + 32; // check if last char is match or mismatch

		targetOfGuide_bit[thr][ti[thr]] = p->splitchar_bit; //
		ti[thr]++;

		d = (guideRNA_bit[guideI[thr]][pos_in_guide] & p->splitchar_bit) != 0 ? d : d - 1; // update distance where d is the number of mismatch

		if (d > -1 && ti[thr] == (len_guide + bulDNA - bD) && // if all the mismatch were used (no more mismatches) and end of the guide
			(max_edits < 0 || (mm - d) + (bulDNA - bD) + (bulRNA - bR) <= max_edits)) // #107 total-edits cap (last char mismatch counted above)
		{ // save results

			save_function_start = omp_get_wtime();

			//saveIndices( p, d, bD, bR, bulType);
			
			saveIndices2(p, d, bD, bR, bulTypeRNA, bulTypeDNA);

			save_function_end = omp_get_wtime();
			save_function += save_function_end - save_function_start;
		}
		d = (guideRNA_bit[guideI[thr]][pos_in_guide] & p->splitchar_bit) != 0 ? d : d + 1; // traceback 
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
	// CRISPRme issue #107: optional 11th arg caps total edits (mismatches + bulges).
	// Omitted or negative => disabled (backward compatible with existing callers).
	if (argc > 11)
		max_edits = atoi(argv[11]);
	// Optional 12th arg: on-the-fly PAM filter (1 = on). Omitted/0 => off. The
	// environment variable CRISPRITZ_PAM_PRUNE (any value) also enables it, so
	// either interface works; the positional arg is what crispritz.py passes for
	// `search --prune-pam`. Kept last so existing 11-arg callers are unaffected.
	if (argc > 12 && atoi(argv[12]) != 0)
		PAM_PRUNE = true;
	// Optional 13th to 15th args, the region hooks (docs/INDEX_FORMAT.md section 13): a
	// skip BED, a keep BED, a leaf-dump file ("" = none; crispritz.py search --skip-bed,
	// --keep-bed, --leaf-dump), each taken from CRISPRITZ_SKIP_BED, CRISPRITZ_KEEP_BED or
	// CRISPRITZ_LEAF_DUMP when its argument is absent or empty. Checked here; the BED is
	// read once the index's contigs are known, before any output is opened.
	const string skipBedPath = regionHookArg(argc, argv, 13, "CRISPRITZ_SKIP_BED");
	const string keepBedPath = regionHookArg(argc, argv, 14, "CRISPRITZ_KEEP_BED");
	leafDumpPath = regionHookArg(argc, argv, 15, "CRISPRITZ_LEAF_DUMP");
	if (!skipBedPath.empty() && !keepBedPath.empty())
	{
		cerr << "ERROR: a search takes a skip BED or a keep BED, not both (skip '" << skipBedPath
			 << "', keep '" << keepBedPath << "')." << endl;
		return 1;
	}
	if (!leafDumpPath.empty() && skipBedPath.empty())
	{
		cerr << "ERROR: the leaf dump ('" << leafDumpPath << "') writes the leaves a skip BED drops, and this search "
			 << "has " << (keepBedPath.empty() ? "no skip BED" : "a keep BED, not a skip BED") << "." << endl;
		return 1;
	}
	regionKeep = !keepBedPath.empty();
	const string regionBedPath = regionKeep ? keepBedPath : skipBedPath;
	regionActive = !regionBedPath.empty();
	{
		const char *sp = getenv("CRISPRITZ_SKIP_PRUNE");
		regionPruneOn = !(sp != nullptr && strcmp(sp, "0") == 0);
	}
	regionPruneActive = regionActive && regionPruneOn;
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

	// ---- On-the-fly PAM filter: one-time init ----
	// Enabled by the 12th positional arg (set above) OR the CRISPRITZ_PAM_PRUNE
	// env var. The search PAM region is `pamRNA` (length L == |pamlimit|). The
	// index stores K flank bases per leaf (K = the build-time pamlimit, encoded
	// in the .bin filename prefix). We support K >= L: a search PAM shorter than
	// the stored flank is N-padded to width K so ONE generous-flank PAM-less index
	// (e.g. an 8-N index) serves ANY PAM of length <= K. The actual padding +
	// searchPamBit construction happen in the length guard below, once K is known
	// (K is read from the first .bin prefix). Here we only flip the mode on.
	// any value but empty or "0": CRISPRITZ_PAM_PRUNE=0 must not switch the filter ON
	const char *prune_env = getenv("CRISPRITZ_PAM_PRUNE");
	if (prune_env != nullptr && prune_env[0] != '\0' && strcmp(prune_env, "0") != 0)
		PAM_PRUNE = true;
	if (PAM_PRUNE)
		cout << "[PAM-FILTER] on-the-fly PAM filter ENABLED; raw search PAM region = "
			 << pamRNA << " (" << pamRNA.size() << " nt)" << endl;
	// ---- END PAM filter init ----

	globalstart = omp_get_wtime(); // start global time
	int numGuide = 0;
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
		if (!pam_at_start)
		{
			reverse(iguide.begin(), iguide.end());
		}

		char *temp_char_guide = (char *)malloc((pamlen - pamlimit + 1) * sizeof(char));
		// guideRNA.push_back((char *)malloc((pamlen - pamlimit) * sizeof(char)));
		guideRNA.push_back(temp_char_guide);
		copy(iguide.begin(), iguide.end(), guideRNA[numGuide]); // save Guide
		guideRNA[numGuide][pamlen - pamlimit] = '\0';
		numGuide++;
		// free(temp_char_guide);
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
			// Only load REGULAR files whose name ends in ".bin" as TST trees. The index
			// directory can also contain provenance sidecars (.pam_build, index_meta.json,
			// .display_label) or subdirectories; handing a non-tree file to loadTST reads
			// garbage (e.g. numNodes=0) and crashes with an out-of-bounds write in
			// deSerialize. Use the same anchored ".bin" suffix test as the PAM-prefix scan
			// below, plus S_ISREG; the >4 guard rejects a bare ".bin" and never admits
			// names like "foo.bin.tmp".
			if (stat(tmp1.c_str(), &file_stats) != 0 || !S_ISREG(file_stats.st_mode) ||
				!(tmp2.size() > 4 && tmp2.compare(tmp2.size() - 4, 4, ".bin") == 0))
				continue;

			file_stats_vec.push_back(file_stats_struct());
			file_stats_vec[file_stats_vec.size() - 1].file_name = dir->d_name;
			file_stats_vec[file_stats_vec.size() - 1].file_dim = (unsigned int)file_stats.st_size;
		}
	}
	closedir(d);

	// A mistyped index path (or an index dir with no .bin trees) must fail loudly
	// rather than silently reporting zero off-targets.
	if (file_stats_vec.empty())
	{
		cerr << "ERROR: no .bin index files found in '" << genome_dir << "'." << endl;
		exit(1);
	}

	// ---- On-the-fly PAM filter: length guard + N-padding to the stored flank ----
	// loadTST reads EXACTLY pamRNA.size() PAM nucleotides per leaf from each .bin,
	// so pamRNA.size() must equal the index's stored flank width K (the build-time
	// pamlimit, encoded in the .bin filename prefix, e.g. "NGG_chr1_0.bin" -> 3,
	// "NNNNNNNN_chr1_0.bin" -> 8). If it did not, the binary stream would be
	// mis-read (the PAM region, the "next" link and every subsequent leaf would
	// shift), silently corrupting results.
	//
	// The "configurable generous flank" idea: build ONE PAM-less index with a wide
	// all-N flank of width K, then serve ANY search PAM of length L <= K by
	// N-PADDING the search PAM up to width K. A real PAM is always adjacent to the
	// protospacer, so the real letters occupy the PROXIMAL positions and the pad N's
	// occupy the DISTAL (outer) positions the index stored but the shorter PAM does
	// not constrain:
	//   * 3' PAM (pam_at_start == false): pamRNA = pam.substr(pamlen - pamlimit,
	//     pamlimit) is ordered [proximal ... distal] (index j increases AWAY from the
	//     protospacer), so the real letters are at low j and we pad on the END:
	//         pamRNA = pamRNA + string(K - L, 'N');
	//   * 5' PAM (pam_at_start == true): pamRNA = pam.substr(0, pamlimit) is ordered
	//     [distal ... proximal] (index j increases TOWARD the protospacer), so the
	//     real letters are at high j and we pad on the FRONT:
	//         pamRNA = string(K - L, 'N') + pamRNA;
	// In both cases the padded pamRNA is laid out in the SAME coordinate frame as the
	// stored flank guideDNA_bit[0..K-1] (loadTST decodes stored position j into
	// guideDNA_bit[j] with the same A=0001/C=0010/G=0100/T=1000 convention used
	// build-side), so searchPamBit[j] compares 1:1 against guideDNA_bit[j] with no
	// reversal, and the pad positions are not tested at all (pamPrune_reject skips
	// them, since a stored 'N' intersects nothing). This is proven
	// equivalent to a dedicated K=L index by test/run_pam_filter_test.sh (an 8-flank
	// NGG search == a dedicated 3-flank NGG index, same result SET).
	//
	// L > K remains a HARD ERROR: the index never stored those extra bases, so a
	// longer search PAM cannot be satisfied.
	// The stored flank width is encoded in the .bin filename prefix. Read it
	// UNCONDITIONALLY: loadTST consumes exactly ceil(pamRNA.size()/2) PAM bytes per
	// leaf, so if that width disagrees with what the index stored, the leaf stream
	// desyncs at the first leaf -- every later PAM, "next" link and finally numNodes
	// is read from a shifted offset. That used to be checked only inside the
	// PAM_PRUNE branch, so a width-K PAM-less index searched with a length-L PAM and
	// the filter OFF walked off the stream instead of being rejected.
	string index_pam_prefix;
	for (size_t f = 0; f < file_stats_vec.size(); ++f)
	{
		const string &fn = file_stats_vec[f].file_name;
		if (fn.size() > 4 && fn.compare(fn.size() - 4, 4, ".bin") == 0)
		{
			size_t us = fn.find('_');
			if (us != string::npos)
			{
				index_pam_prefix = fn.substr(0, us);
				break;
			}
		}
	}
	// So every .bin must be named for one PAM. Files named for another are not part of
	// the index: a macOS "._" (AppleDouble) file copied along with an index is one, and
	// listed after the index's own files it used to be read as a tree (its first int as
	// a leaf count of 118,883,584). Refused whatever order the directory lists them in,
	// as crispritz.py index-info refuses such a directory.
	{
		vector<string> pams;
		for (size_t f = 0; f < file_stats_vec.size(); ++f)
		{
			const string &fn = file_stats_vec[f].file_name;
			string named = fn.substr(0, fn.find('_'));
			if (find(pams.begin(), pams.end(), named) == pams.end())
				pams.push_back(named);
		}
		if (pams.size() > 1)
		{
			sort(pams.begin(), pams.end());
			cerr << "ERROR: the .bin files in '" << genome_dir << "' are named for more than one PAM (";
			for (size_t k = 0; k < pams.size(); ++k)
				cerr << (k ? ", '" : "'") << pams[k] << "'";
			cerr << "): an index directory holds the .bin files of one index only (a macOS \"._\" "
				 << "file copied along with an index is not one of them); remove the others." << endl;
			return 1;
		}
	}
	if (!PAM_PRUNE && !index_pam_prefix.empty()
		&& index_pam_prefix.size() != pamRNA.size())
	{
		cerr << "ERROR: this index stores a " << index_pam_prefix.size()
			 << "-nt PAM/flank region ('" << index_pam_prefix << "') but the search "
			 << "PAM region '" << pamRNA << "' is " << pamRNA.size() << " nt. Reading "
			 << "it at the wrong width would mis-parse every leaf of the index, so the "
			 << "search is refused." << endl;
		if (index_pam_prefix.find_first_not_of('N') == string::npos)
			cerr << "       This is a PAM-less index: pass --prune-pam to serve a "
				 << "shorter PAM from it, or use an index built for this PAM." << endl;
		else
			cerr << "       Use an index built for this PAM." << endl;
		return 1;
	}
	if (PAM_PRUNE)
	{
		if (!index_pam_prefix.empty())
		{
			size_t K = index_pam_prefix.size(); // stored flank width in the index
			size_t L = pamRNA.size();            // requested search PAM length
			if (L > K)
			{
				cerr << "ERROR: --prune-pam/CRISPRITZ_PAM_PRUNE search PAM is LONGER "
					 << "than the index's stored flank. Index flank is '"
					 << index_pam_prefix << "' (" << K << " nt) but search PAM region "
					 << "is '" << pamRNA << "' (" << L << " nt). A longer PAM cannot "
					 << "be served because the index never stored those extra flank "
					 << "bases; rebuild the index with an all-N PAM of length >= "
					 << L << " (recommended generous default: 8)." << endl;
				return 1;
			}
			if (L < K)
			{
				// N-pad the search PAM up to the stored flank width K. Pad DISTAL:
				// end for a 3' PAM, front for a 5' PAM (see rationale above).
				if (pam_at_start)
				{
					// The emit path trims prune_pad from the END of the reported
					// crRNA/DNA (see the `pamRNA.size() - prune_pad` sites). For a 5'
					// PAM the pad goes on the FRONT, so that trim removes real bases
					// and leaves the pad in the reported PAM region -- loci still
					// match a dedicated index, which is all the equivalence test
					// checks, but every emitted row carries N's where the PAM should
					// be and a downstream PAM filter then drops all of them. Refuse
					// until the 5' emit/coordinate reconciliation is implemented.
					cerr << "ERROR: --prune-pam cannot yet serve a SHORTER PAM from a "
						 << "wider flank for a 5' PAM. Index flank is '"
						 << index_pam_prefix << "' (" << K << " nt), search PAM '"
						 << pamRNA << "' (" << L << " nt), PAM at start. The emitted "
						 << "sequence is trimmed at the end only, so a front pad would "
						 << "be reported inside the PAM region. Use an index whose "
						 << "flank width equals the PAM length (" << L << ") for a 5' "
						 << "PAM." << endl;
					return 1;
				}
				string pad(K - L, 'N');
				prune_pad = (int)(K - L); // minus-strand (3') coord reconciliation to a dedicated L index
				pamRNA = pamRNA + pad;   // 3': real letters proximal (low j)
				cout << "[PAM-FILTER] N-padded search PAM to the index flank width K="
					 << K << " (orientation: " << (pam_at_start ? "5'/front" : "3'/end")
					 << "-pad): effective PAM region = " << pamRNA << endl;
			}
		}
		// Build the IUPAC-expanded bitsets now that pamRNA is at the stored width K.
		searchPamBit = pamPrune_bit_conversion(pamRNA);
		cout << "[PAM-FILTER] effective search PAM region = " << pamRNA << " ("
			 << searchPamBit.size() << " nt)" << endl;
	}
	// ---- END PAM filter length guard + N-padding ----

	// The build PAM, for recognising the format-1 index-0 artefact in loadTST.
	buildPamBit = pamPrune_bit_conversion(index_pam_prefix);
	for (size_t j = 0; j < buildPamBit.size(); j++)
		if (!buildPamBit[j].all())
			buildPamConstrained = true;

	// Read every .bin's header before searching (docs/INDEX_FORMAT.md): refuse a
	// format this engine does not know, a malformed or cut-short format-2 file, and a
	// format-2 file built for another PAM width or PAM side than the search reads it
	// with -- before any thread starts. A format-1 file (CRISPRitz <= 2.9.1) has no
	// header and is searched as 2.9.1 searched it, but for the index-0 fix (IDX-14).
	// A search guide of another length than the index's is searched (as in 2.9.1),
	// with a note on stdout saying what that changes (section 5.1).
	{
		int nFormat[3] = {0, 0, 0};
		int nQuad = 0;
		bool mbNoted = false;
		vector<int> guideLenNoted;
		for (size_t f = 0; f < file_stats_vec.size(); ++f)
		{
			const string &fn = file_stats_vec[f].file_name;
			string path = string(genome_dir) + "/" + fn;
			ifstream hf(path, ios::in | ios::binary);
			IndexHeader h;
			string err;
			if (!hf)
				err = "cannot be opened";
			if (!err.empty() || !readIndexHeader(hf, h, err))
			{
				cerr << "ERROR: '" << path << "' " << err << "." << endl;
				return 1;
			}
			nFormat[h.format]++;
			if (h.format == 2 && (h.flags & CZI_FLAG_FOUR_ALLELE))
				nQuad++;
			// the body opens with numLeaves and the guide length the index was built
			// for (a file too short for them is refused by loadTST, as in 2.9.1)
			int body[2];
			bool haveBody = false;
			if (hf.read((char *)body, sizeof(body)))
			{
				haveBody = true;
				if (h.format == 2 && body[1] != h.guideLen)
				{
					cerr << "ERROR: '" << path << "' has a malformed index: its header gives guide length "
						 << h.guideLen << ", its body " << body[1] << "." << endl;
					return 1;
				}
				int gi = body[1], gs = pamlen - pamlimit, imb = h.format == 2 ? h.mb : max_bulges;
				if (gi != gs && find(guideLenNoted.begin(), guideLenNoted.end(), gi) == guideLenNoted.end())
				{
					// stdout, not stderr: CRISPRme treats any stderr output as a failure
					guideLenNoted.push_back(gi);
					cout << "NOTE: the index was built for " << gi << "-nt guides and the search guides are "
						 << gs << " nt; ";
					if (gs > gi)
						cout << "it holds " << gi + imb << " bases a site, so an alignment is found only if its "
							 << "DNA bulges minus RNA bulges are at most " << gi + imb - gs << " (an index built for "
							 << gs << "-nt guides with max_bulges " << imb << " allows " << imb << ").";
					else
						cout << "it holds a site only if the site's " << gi << "-nt guide region lies on the "
							 << "contig with no N, so a target within " << gi - gs << " nt of a contig end or "
							 << "an N run on its side away from the PAM is not found (an index built for " << gs
							 << "-nt guides finds it).";
					cout << " See docs/INDEX_FORMAT.md, section 5.1." << endl;
				}
			}
			string prefix = fn.substr(0, fn.find('_'));
			if (h.format == 2)
			{
				if (h.pam != prefix || (size_t)h.pamLen != pamRNA.size())
				{
					cerr << "ERROR: '" << path << "' holds a " << h.pamLen << "-nt PAM/flank region ('"
						 << h.pam << "') but is named for '" << prefix << "' and searched at "
						 << pamRNA.size() << " nt; refusing to read it." << endl;
					return 1;
				}
				if (((h.flags & CZI_FLAG_PAM_AT_START) != 0) != pam_at_start)
				{
					cerr << "ERROR: '" << path << "' was built for a PAM at the "
						 << ((h.flags & CZI_FLAG_PAM_AT_START) ? "5'" : "3'") << " end, but the search PAM is at the "
						 << (pam_at_start ? "5'" : "3'") << " end (the sign of the PAM file's length)." << endl;
					return 1;
				}
				if (h.mb != max_bulges && !mbNoted)
				{
					// stdout, not stderr: CRISPRme treats any stderr output as a failure
					cout << "NOTE: the index was built with max_bulges=" << h.mb << " (command line: "
						 << max_bulges << "); its coordinates are decoded with " << h.mb << "." << endl;
					mbNoted = true;
				}
			}
			// A leaf record takes at least 4 + ceil(K/2) + 1 bytes (the window start, the
			// flank nibbles, the marker or flags byte), so a leaf count the file cannot
			// hold is damage (a format-1 file has no end marker to show it): refused here,
			// where loadTST would allocate that many leaves and fail in std::bad_alloc.
			if (haveBody && body[0] > 0)
			{
				long long width = h.format == 2 ? h.pamLen : (long long)prefix.size();
				long long need = (h.format == 2 ? h.headerLen : 0) + 12 + (long long)body[0] * (4 + (width + 1) / 2 + 1);
				hf.clear();
				hf.seekg(0, ios::end);
				long long size = (long long)hf.tellg();
				if (need > size)
				{
					cerr << "ERROR: '" << path << "' is damaged: its " << body[0] << " leaves need at least "
						 << need << " bytes, the file holds " << size << "." << endl;
					return 1;
				}
			}
		}
		cout << "INDEX FORMAT:";
		if (nFormat[1])
			cout << " 1 (CRISPRitz <= 2.9.1) x" << nFormat[1];
		if (nFormat[2])
			cout << " 2 x" << nFormat[2];
		if (nQuad)
			cout << " (" << nQuad << " holding Q, the four-allele symbol)";
		cout << endl;
	}

	// Region hooks: read the BED against the index's contigs (those its .bin files are
	// named for), and create the leaf dump, before any output is opened.
	if (regionActive)
	{
		set<string> contigs;
		for (size_t f = 0; f < file_stats_vec.size(); ++f)
			contigs.insert(binContig(file_stats_vec[f].file_name));
		string err;
		long long nIntervals = 0;
		if (!loadRegionBed(regionBedPath, contigs, err, nIntervals))
		{
			cerr << "ERROR: the " << (regionKeep ? "keep" : "skip") << " BED '" << regionBedPath << "' " << err
				 << "; the search is refused (docs/INDEX_FORMAT.md section 13)." << endl;
			return 1;
		}
		if (!leafDumpPath.empty())
		{
			ofstream df(leafDumpPath, ios::out | ios::trunc);
			if (!df)
			{
				cerr << "ERROR: the leaf dump '" << leafDumpPath << "' cannot be written." << endl;
				return 1;
			}
		}
		size_t nMerged = 0;
		for (const auto &kv : regionIv)
			nMerged += kv.second.size();
		cout << (regionKeep ? "KEEP BED " : "SKIP BED ") << regionBedPath << ": " << nIntervals << " interval(s), "
			 << nMerged << " once merged, on " << regionIv.size() << " of the index's " << contigs.size()
			 << " contig(s); " << (regionKeep ? "only the sites whose core meets one are searched"
											  : "the sites whose core meets one are skipped")
			 << " (subtree pruning " << (regionPruneOn ? "on" : "off") << ")" << endl;
	}

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
		create_profile = true;
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
	profiling_rna_mm.resize(pamlen - pamlimit + (mm + 1) * max_bulges, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));
	profiling_rna.resize(pamlen - pamlimit + (mm + 1) * max_bulges, vector<vector<int>>(numGuide, vector<int>(num_thr, 0)));



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
	idxFacts.resize(num_thr);
	if (regionActive)
	{
		regionTree.assign(num_thr, 0);
		nodeLive.resize(num_thr);
		nRegionSkipped.assign(num_thr, 0);
		nRegionKept.assign(num_thr, 0);
		leafDump.resize(num_thr);
	}
	albero_glb.resize(num_thr);
	ti.resize(num_thr);
	gi.resize(num_thr);
	inGuide_bit.resize(num_thr);
	targetOfGuide_bit.resize(num_thr);
	len_guide = pamlen - pamlimit;
	inGuide.resize(num_thr);
	targetOfGuide.resize(num_thr);

	//Matrix for joint DNA+RNA profiling
	int bins_both = len_guide + bulDNA + (mm + 1) * (max_bulges + 1) * (max_bulges + 1);
	profiling_both_mm.resize(
		bins_both,
		vector<vector<int>>(numGuide, vector<int>(num_thr, 0))
	);

	int numNodes;
	int numLeaves;

	char c_inGuide[100];
	char c_targetOfGuide[100];
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
			save_function = 0;
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

			double guidecut = 200; //50;
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


					nearsearch(guideRNA[i], &albero_glb[thr][0], 0, mm2, bulDNA2, bulRNA2, true, 0, 0, thr);
					
					
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

							fileResults << bulgeType[thr][j] << "\t" << vecInGuide[thr][j] << "\t" << vecTargetOfGuide[thr][j] << "\t" << chrName_glb[thr][j] << "\t" << indices[thr][j] << "\t" << cluster_position[thr][j] << "\t" << directions[thr][j] << "\t" << mismatches[thr][j] << "\t" << bulgeSize[thr][j] << "\t" << mismatches[thr][j] + bulgeSize[thr][j] << "\n";
							// cout << bulgeType[thr][j] << "\t" << vecInGuide[thr][j] << "\t" << vecTargetOfGuide[thr][j] << "\t" << chrName_glb[thr][j] << "\t" << indices[thr][j] << "\t" << cluster_position[thr][j] << "\t" << directions[thr][j] << "\t" << mismatches[thr][j] << "\t" << bulgeSize[thr][j] << "\t" << mismatches[thr][j] + bulgeSize[thr][j] << "\n";
						}
					}

					bulgeType[thr].clear();		   //bulgeType.resize(num_thr);
					vecInGuide[thr].clear();	   //vecInGuide.resize(num_thr);
					vecTargetOfGuide[thr].clear(); //vecTargetOfGuide.resize(num_thr);
					chrName_glb[thr].clear();	   //chrName_glb.resize(num_thr);
					indices[thr].clear();		   //indices.resize(num_thr);
					cluster_position[thr].clear();
					directions[thr].clear(); //directions.resize(num_thr);
					mismatches[thr].clear(); //mismatches.resize(num_thr);
					bulgeSize[thr].clear();	 //bulgeSize.resize(num_thr);
				}
			} //uncomment/comment to activate/deactivate saving every 50 guides

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
		cout << "ENTERED PROFILE" << "\n";

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
			for (int j = 1; j <= max_bulges; j++)
			{
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
			for (int j = 1; j <= max_bulges; j++)
			{
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
			saveProfileGuide(
				guideRNA_s[i], i,
				mm, max_bulges, len_guide, bulDNA,
				profiling,               // base (pos & mm-only summary)
				ext_profiling,           // extended nuc×pos per mm
				profiling_dna, profiling_dna_mm,
				profiling_rna, profiling_rna_mm,
				profiling_both_mm,       // <-- NEW joint bins (mm, bulgeDNA, bulgeRNA)
				fileprofiling, file_ext_profiling,
				file_profiling_dna, file_profiling_rna, file_profiling_complete,
				num_thr, (int)pamRNA.size(), pam_at_start,
				ext_profiling_dna, ext_profiling_rna
			);
		}
		fileprofiling.close();
		file_ext_profiling.close();
		file_profiling_dna.close();
		file_profiling_rna.close();
		file_profiling_complete.close();
	}

	fileResults.close();

	// Remove duplicated lines in targets.txt file
	removeDuplicatesInPlace(name_result + ".targets.txt");

	// Region hooks: the leaves the BED dropped and kept (CRISPRme requires this line of
	// a search it gave a BED), and the leaf dump, sorted, one line a leaf.
	if (regionActive)
	{
		long long nSkipped = 0, nKept = 0;
		for (size_t t = 0; t < nRegionSkipped.size(); t++)
		{
			nSkipped += nRegionSkipped[t];
			nKept += nRegionKept[t];
		}
		cout << "SKIPPED LEAVES\t" << nSkipped << "\tKEPT LEAVES\t" << nKept << "\n";
		if (!leafDumpPath.empty())
		{
			vector<LeafDumpLine> all;
			for (auto &b : leafDump)
			{
				all.insert(all.end(), b.begin(), b.end());
				vector<LeafDumpLine>().swap(b);
			}
			sort(all.begin(), all.end()); // one line a dropped leaf (a 2.9.1 PAM-less index holds index 0 twice)
			ofstream df(leafDumpPath, ios::out | ios::trunc);
			for (const LeafDumpLine &d : all)
				df << d.contig << "\t" << d.ws << "\t" << d.letters << "\t" << d.side << "\n";
			df.close();
			if (!df)
			{
				cerr << "ERROR: the leaf dump '" << leafDumpPath << "' could not be written in full." << endl;
				return 1;
			}
			cout << "LEAF DUMP\t" << all.size() << "\n";
		}
	}
	globalend = omp_get_wtime(); // end global time
	cout << "-----------------------"
		 << "\n";
	cout << "TOTAL TIME:\t" << globalend - globalstart << "\n";

	return 0;
}

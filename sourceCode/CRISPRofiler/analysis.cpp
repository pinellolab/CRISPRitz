#include "include/crispritz.h"

using namespace std;

vector<int> pamindices, pamindicesreverse;
extern vector<bitset<4>> genomebit;
extern string genome, pam;
extern int pamlen, pamlimit;
extern bool pamdirection;

void analysis()
{
	//variable reset
	pamindices.clear();
	pamindicesreverse.clear();

	//start pam searching
	// searchPam();
	searchPAMonGenome(pam, pamlen, genome, pamlimit, pamdirection, 0, 0);
	// vector<int> searchPAMonGenome(string pam_sequence, int pam_len, string genome_sequence, int pam_limit, bool pam_at_start, int max_bulges, int max_mismatches)

	//start guides searching, executed number of guides times
	guide_searching();

	//clear the genome string for the next genome analysis
	genome.clear();
	genomebit.clear();
}

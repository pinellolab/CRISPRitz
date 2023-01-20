#include "include/crispritz.h"

using namespace std;

vector<int> pamindices, pamindicesreverse;
extern vector<bitset<4>> genomebit;
extern string genome, pam;
extern int pamlen, pamlimit, guidelen;
extern bool pam_at_start;

void analysis()
{
	// variable reset
	pamindices.clear();
	pamindicesreverse.clear();

	// start pam searching
	searchPAMonGenome(pam, guidelen, genome, pamlimit, pam_at_start, 0, 0);

	// start guides searching, executed number of guides times
	guide_searching();

	// clear the genome string for the next genome analysis
	genome.clear();
	genomebit.clear();
}

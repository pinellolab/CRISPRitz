#include "include/crispritz.h"

using namespace std;

vector<int> pamindices, pamindicesreverse;
extern vector<bitset<4>> genomebit;
extern string genome;

void analysis()
{
	//variable reset
	pamindices.clear();
	pamindicesreverse.clear();

	//start pam searching
	searchPam();

	//start guides searching, executed number of guides times
	guide_searching();

	//clear the genome string for the next genome analysis
	genome.clear();
	genomebit.clear();
}
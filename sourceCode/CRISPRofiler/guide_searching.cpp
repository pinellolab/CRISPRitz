#include "include/crispritz.h"

using namespace std;

int tmp_mismatch;
string totalbuf;
extern ofstream results;
extern int currentmissmatch, pamlimit, guidelen, guidecount, totalguides, threads, inputmissmatch;
extern vector<int> missmatchthrestotal, pamindices, pamindicesreverse;
extern vector<string> guides;
extern string chrnames, genome;
extern char nowriting, noprofile;
extern bool pam_at_start;

// profile
extern vector<vector<int>> guideprofiling;
extern vector<vector<vector<vector<int>>>> matrixprofiling;

// genome
extern vector<bitset<4>> genomebit;
extern vector<vector<bitset<4>>> guidesbit;
extern vector<vector<bitset<4>>> reverseguidesbit;

string reversetarget(string targetfound) // reverse del target per poterlo confrontare con la guida
{
	// #pragma omp parallel for
	string complement = targetfound;
	vector<int> nt_to_lower;

	// cout << "target found prima " << targetfound << endl;

	for (int nt = 0; nt < targetfound.length(); nt++)
	{
		if (islower(targetfound[nt]))
			nt_to_lower.push_back(nt);
		complement[nt] = toupper(complement[nt]);
		targetfound[nt] = toupper(targetfound[nt]);
	}

	// cout << "target found dopo " << targetfound << " complement dopo " << complement << endl;

	for (int j = 0; j < targetfound.length(); j++)
	{
		if (targetfound[j] == 'A')
		{
			complement[j] = 'T';
		}
		else if (targetfound[j] == 'T')
		{
			complement[j] = 'A';
		}
		else if (targetfound[j] == 'G')
		{
			complement[j] = 'C';
		}
		else if (targetfound[j] == 'C')
		{
			complement[j] = 'G';
		}
		else if (targetfound[j] == 'R')
		{
			complement[j] = 'Y';
		}
		else if (targetfound[j] == 'Y')
		{
			complement[j] = 'R';
		}
		else if (targetfound[j] == 'S')
		{
			complement[j] = 'S';
		}
		else if (targetfound[j] == 'W')
		{
			complement[j] = 'W';
		}
		else if (targetfound[j] == 'M')
		{
			complement[j] = 'K';
		}
		else if (targetfound[j] == 'K')
		{
			complement[j] = 'M';
		}
		else if (targetfound[j] == 'H')
		{
			complement[j] = 'D';
		}
		else if (targetfound[j] == 'D')
		{
			complement[j] = 'H';
		}
		else if (targetfound[j] == 'B')
		{
			complement[j] = 'V';
		}
		else if (targetfound[j] == 'V')
		{
			complement[j] = 'B';
		}
		else if (targetfound[j] == 'N')
		{
			complement[j] = 'N';
		}
	}

	for (int nt = 0; nt < nt_to_lower.size(); nt++)
		complement[nt_to_lower[nt]] = tolower(complement[nt_to_lower[nt]]);

	reverse(complement.begin(), complement.end());

	// cout << "target found dopo dopo " << targetfound << " complement dopo dopo " << complement << endl;

	return complement;
}

vector<bitset<4>> string_to_bit_converter(string sequence)
{
	vector<bitset<4>> sequence_bit;
	sequence_bit.resize(sequence.length());

	for (int i = 0; i < sequence.length(); ++i)
	{
		if (sequence[i] == 'A')
		{
			sequence_bit[i] = bitset<4>(string("0001"));
		}
		else if (sequence[i] == 'C')
		{
			sequence_bit[i] = bitset<4>(string("0010"));
		}
		else if (sequence[i] == 'G')
		{
			sequence_bit[i] = bitset<4>(string("0100"));
		}
		else if (sequence[i] == 'T')
		{
			sequence_bit[i] = bitset<4>(string("1000"));
		}
		else if (sequence[i] == 'N')
		{
			sequence_bit[i] = bitset<4>(string("1111"));
		}
		else if (sequence[i] == 'R')
		{
			sequence_bit[i] = bitset<4>(string("0101"));
		}
		else if (sequence[i] == 'Y')
		{
			sequence_bit[i] = bitset<4>(string("1010"));
		}
		else if (sequence[i] == 'S')
		{
			sequence_bit[i] = bitset<4>(string("0110"));
		}
		else if (sequence[i] == 'W')
		{
			sequence_bit[i] = bitset<4>(string("1001"));
		}
		else if (sequence[i] == 'K')
		{
			sequence_bit[i] = bitset<4>(string("1100"));
		}
		else if (sequence[i] == 'M')
		{
			sequence_bit[i] = bitset<4>(string("0011"));
		}
		else if (sequence[i] == 'B')
		{
			sequence_bit[i] = bitset<4>(string("1110"));
		}
		else if (sequence[i] == 'D')
		{
			sequence_bit[i] = bitset<4>(string("1101"));
		}
		else if (sequence[i] == 'H')
		{
			sequence_bit[i] = bitset<4>(string("1011"));
		}
		else if (sequence[i] == 'V')
		{
			sequence_bit[i] = bitset<4>(string("0111"));
		}
	}
	return sequence_bit;
}

string missmatching(string targetfound, int index, int guidafound, int reverse) // conteggio dei mismatches e profilazione se richiesta dall'utente
{
	// assegnazione valori ai contatori temporanei
	vector<int> tmp_mismatch_position;
	vector<vector<float>> tmp_bp_mismatch;
	vector<bitset<4>> sequence_bit;
	tmp_mismatch = 0;
	tmp_mismatch_position.resize(targetfound.length());
	tmp_bp_mismatch.resize(targetfound.length());

	// cout << "entro in missmatching con seq " << targetfound << " e guida " << guides[guidafound] << " in pos " << index << " reverse è " << reverse << endl;
	sequence_bit = string_to_bit_converter(targetfound);

	for (int j = 0; j < targetfound.length(); j++)
	{
		tmp_bp_mismatch[j].resize(4); // contatore per motiflogo
		// check if mm in bp genome/guide
		if ((sequence_bit[j] & guidesbit[guidafound][j]) == 0)
		{
			if (targetfound[j] == 'A')
				tmp_bp_mismatch[j][0]++;
			else if (targetfound[j] == 'T')
				tmp_bp_mismatch[j][3]++;
			else if (targetfound[j] == 'G')
				tmp_bp_mismatch[j][2]++;
			else if (targetfound[j] == 'C')
				tmp_bp_mismatch[j][1]++;
			else if (targetfound[j] == 'R')
			{
				tmp_bp_mismatch[j][0] = tmp_bp_mismatch[j][0] + 0.5;
				tmp_bp_mismatch[j][2] = tmp_bp_mismatch[j][2] + 0.5;
			}
			else if (targetfound[j] == 'Y')
			{
				tmp_bp_mismatch[j][1] = tmp_bp_mismatch[j][1] + 0.5;
				tmp_bp_mismatch[j][3] = tmp_bp_mismatch[j][3] + 0.5;
			}
			else if (targetfound[j] == 'S')
			{
				tmp_bp_mismatch[j][2] = tmp_bp_mismatch[j][2] + 0.5;
				tmp_bp_mismatch[j][1] = tmp_bp_mismatch[j][1] + 0.5;
			}
			else if (targetfound[j] == 'W')
			{
				tmp_bp_mismatch[j][0] = tmp_bp_mismatch[j][0] + 0.5;
				tmp_bp_mismatch[j][3] = tmp_bp_mismatch[j][3] + 0.5;
			}
			else if (targetfound[j] == 'K')
			{
				tmp_bp_mismatch[j][2] = tmp_bp_mismatch[j][2] + 0.5;
				tmp_bp_mismatch[j][3] = tmp_bp_mismatch[j][3] + 0.5;
			}
			else if (targetfound[j] == 'M')
			{
				tmp_bp_mismatch[j][0] = tmp_bp_mismatch[j][0] + 0.5;
				tmp_bp_mismatch[j][1] = tmp_bp_mismatch[j][1] + 0.5;
			}
			else if (targetfound[j] == 'B')
			{
				tmp_bp_mismatch[j][1] = tmp_bp_mismatch[j][1] + (1 / 3);
				tmp_bp_mismatch[j][2] = tmp_bp_mismatch[j][2] + (1 / 3);
				tmp_bp_mismatch[j][3] = tmp_bp_mismatch[j][3] + (1 / 3);
			}
			else if (targetfound[j] == 'D')
			{
				tmp_bp_mismatch[j][0] = tmp_bp_mismatch[j][0] + (1 / 3);
				tmp_bp_mismatch[j][2] = tmp_bp_mismatch[j][2] + (1 / 3);
				tmp_bp_mismatch[j][3] = tmp_bp_mismatch[j][3] + (1 / 3);
			}
			else if (targetfound[j] == 'H')
			{
				tmp_bp_mismatch[j][0] = tmp_bp_mismatch[j][0] + (1 / 3);
				tmp_bp_mismatch[j][1] = tmp_bp_mismatch[j][1] + (1 / 3);
				tmp_bp_mismatch[j][3] = tmp_bp_mismatch[j][3] + (1 / 3);
			}
			else if (targetfound[j] == 'V')
			{
				tmp_bp_mismatch[j][0] = tmp_bp_mismatch[j][0] + (1 / 3);
				tmp_bp_mismatch[j][1] = tmp_bp_mismatch[j][1] + (1 / 3);
				tmp_bp_mismatch[j][2] = tmp_bp_mismatch[j][2] + (1 / 3);
			}

			targetfound[j] = tolower(targetfound[j]); // lettera minuscola
			tmp_mismatch++;							  // conto MM totali della guida
			guideprofiling[guidafound][j]++;		  // segno un MM in una base della guida nel profilo
			tmp_mismatch_position[j] = 5;			  // segno posizione MM
		}
	}

	// cout << "faccio mismatching " << endl;

	if (tmp_mismatch > 0)
	{
		// aggiorno il profili della guida con OFF-TARGET E MM TOTALI
		guideprofiling[guidafound][guidelen + 1]++;
		guideprofiling[guidafound][guidelen + 2] += tmp_mismatch;
	}
	else
	{
		// aggiorno profilo con ON-TARGET
		guideprofiling[guidafound][guidelen]++;
	}

	// cout << "faccio profiling 1 " << endl;

	// aggiorno profilo con nuovi MM nella colonna con threshold
	guideprofiling[guidafound][guidelen + 3 + tmp_mismatch]++;

	// cout << "faccio profiling 2 " << endl;

	// aggiorno extendend_profile con dati della guida
	for (int gg = 0; gg < targetfound.length(); gg++)
	{
		if (tmp_mismatch_position[gg] == 5)
		{
			// cout << "matrixprof " << matrixprofiling[guidafound][tmp_mismatch][gg][4] << endl;
			matrixprofiling[guidafound][tmp_mismatch][gg][4]++;
		}
		for (int dd = 0; dd < 4; dd++)
		{
			// cout << "matrixprof " << matrixprofiling[guidafound][tmp_mismatch][gg][dd] << endl;
			matrixprofiling[guidafound][tmp_mismatch][gg][dd] += tmp_bp_mismatch[gg][dd];
		}
	}

	// cout << "faccio matrix profiling" << endl;

	return targetfound;
}

// execute guide search
void guide_searching()
{
	vector<int> respos, resneg, respos_private, resneg_private;			// inizializzazione variabili per store indici pam
	vector<int> guidepos, guideneg, guidepos_private, guideneg_private; // inizializzaione variabili per store indici guide

	// count the number of indices found in each PAM direction
	int pampossize = pamindices.size();
	int pamnegsize = pamindicesreverse.size();

	// cout << "entro nel search guide" << endl;

	// for (int guide = 0; guide < guides.size(); guide++)
	// {
	// 	for (int gg = 0; gg < guides[guide].length(); gg++)
	// 	{
	// 		cout << "guida is " << guides[guide] << endl;
	// 		cout << "bit is " << guidesbit[guide][gg] << endl;
	// 		cout << "reverse bit is " << reverseguidesbit[guide][gg] << endl;
	// 	}
	// 	cout << "fine guida" << endl;
	// }

	// parallel region to find all targets on a selected chromosome
#pragma omp parallel num_threads(threads) private(respos_private, guidepos_private, resneg_private, guideneg_private)
	{
#pragma omp for simd schedule(static) nowait
		for (int i = 0; i < pampossize; i++) // ciclo sugli indici trovati precendetemente, PAM
		{
			for (int guidecount = 0; guidecount < guides.size(); guidecount++) // ciclo le guide, su tutti gli indici
			{
				int currentmissmatch = 0; // missmatch di ogni guida

				for (int j = 0; j < guides[guidecount].length() & (currentmissmatch <= inputmissmatch); j++) // confronto lettera-lettera, guida su genoma
				{
					if ((genomebit[j + pamindices[i] - (guides[guidecount].length() - pamlimit)] & guidesbit[guidecount][j]) == 0)
					{
						currentmissmatch++;
						if (genomebit[j + pamindices[i]] == 0)
							currentmissmatch += 100;
					}
				}
				if (currentmissmatch <= inputmissmatch) // se rimango nella soglia, salvo dati, guida e pam che hanno generato target
				{
					respos_private.push_back(pamindices[i] - (guides[guidecount].length() - pamlimit));
					guidepos_private.push_back(guidecount);
				}
			}
		}

#pragma omp for simd schedule(static) nowait
		for (int i = 0; i < pamnegsize; i++)
		{
			for (int guidecount = 0; guidecount < guides.size(); guidecount++)
			{
				int currentmissmatch = 0;

				for (int j = 0; j < guides[guidecount].length() & (currentmissmatch <= inputmissmatch); j++)
				{
					if ((genomebit[j + pamindicesreverse[i]] & reverseguidesbit[guidecount][j]) == 0)
					{
						currentmissmatch++;
						if (genomebit[j + pamindicesreverse[i]] == 0)
							currentmissmatch += 100;
					}
				}
				if (currentmissmatch <= inputmissmatch)
				{
					resneg_private.push_back(pamindicesreverse[i]);
					guideneg_private.push_back(guidecount);
				}
			}
		}
#pragma omp critical
		{
			// copy all threads values in a single value usable to write results
			respos.insert(respos.end(), respos_private.begin(), respos_private.end());
			guidepos.insert(guidepos.end(), guidepos_private.begin(), guidepos_private.end());

			resneg.insert(resneg.end(), resneg_private.begin(), resneg_private.end());
			guideneg.insert(guideneg.end(), guideneg_private.begin(), guideneg_private.end());
		}
	}

	// cout << "found results pos " << respos.size() << " found results neg" << resneg.size() << endl;

	// if (!pam_at_start) // direzione pam
	// {
	// check delle posizioni con risultati positivi
	for (int i = 0; i < respos.size(); i++)
	{
		totalbuf += "X"; // bulge type (X stands for no bulge)
		totalbuf += "\t";
		if (!pam_at_start)
			totalbuf += guides[guidepos[i]]; // guida
		else
			totalbuf += reversetarget(guides[guidepos[i]]); // guida
		totalbuf += "\t";
		if (!pam_at_start)
			totalbuf += missmatching(genome.substr(respos[i], guides[guidepos[i]].length()), respos[i], guidepos[i], 0); // calcolo della stringa e aggiornamento contatori profili
		else
			totalbuf += reversetarget(missmatching(genome.substr(respos[i], guides[guidepos[i]].length()), respos[i], guidepos[i], 0)); // calcolo della stringa e aggiornamento contatori profili
		totalbuf += "\t";
		totalbuf += chrnames; // nome chr
		totalbuf += "\t";
		totalbuf += to_string(respos[i]); // posizione del TARGET nel genoma
		totalbuf += "\t";
		totalbuf += to_string(respos[i]); // Cluster Position del TARGET nel genoma
		totalbuf += "\t";
		if (!pam_at_start)
			totalbuf += "+"; // strand
		else
			totalbuf += "-"; // strand
		totalbuf += "\t";
		totalbuf += to_string(tmp_mismatch); // mismatch per target
		totalbuf += "\t";
		totalbuf += "0"; // num bulge, sempre uguale 0
		totalbuf += "\t";
		totalbuf += to_string(tmp_mismatch); // total per target
		totalbuf += "\n";
	}

	// check delle posizioni con risultati negativi
	for (int i = 0; i < resneg.size(); i++)
	{
		string target = reversetarget(genome.substr(resneg[i], guides[guideneg[i]].length()));
		totalbuf += "X"; // bulge type (X stands for no bulge)
		totalbuf += "\t";
		if (!pam_at_start)
			totalbuf += guides[guideneg[i]]; // guida
		else
			totalbuf += reversetarget(guides[guideneg[i]]); // guida
		totalbuf += "\t";
		if (!pam_at_start)
			totalbuf += missmatching(target, resneg[i], guideneg[i], 1);
		else
			totalbuf += reversetarget(missmatching(target, resneg[i], guideneg[i], 1));
		totalbuf += "\t";
		totalbuf += chrnames;
		totalbuf += "\t";
		totalbuf += to_string(resneg[i]);
		totalbuf += "\t";
		totalbuf += to_string(resneg[i]); // Cluster Position
		totalbuf += "\t";
		if (!pam_at_start)
			totalbuf += "-";
		else
			totalbuf += "+";
		totalbuf += "\t";
		totalbuf += to_string(tmp_mismatch);
		totalbuf += "\t";
		totalbuf += "0";
		totalbuf += "\t";
		totalbuf += to_string(tmp_mismatch); // total per target
		totalbuf += "\n";
	}
	// }
	// else
	// {
	// 	for (int i = 0; i < respos.size(); i++)
	// 	{
	// 		// string target = reversetarget(genome.substr(respos[i], guides[guidepos[i]].length()));
	// 		totalbuf += "X"; // bulge type (X stands for no bulge)
	// 		totalbuf += "\t";
	// 		totalbuf += guides[guidepos[i]]; // guida
	// 		totalbuf += "\t";
	// 		totalbuf += missmatching(genome.substr(respos[i], guides[guidepos[i]].length()), respos[i], guidepos[i], 0); // calcolo della stringa e aggiornamento contatori profili
	// 		totalbuf += "\t";
	// 		totalbuf += chrnames; // nome chr
	// 		totalbuf += "\t";
	// 		totalbuf += to_string(respos[i]); // posizione del TARGET nel genoma
	// 		totalbuf += "\t";
	// 		totalbuf += to_string(respos[i]); // Cluster Position del TARGET nel genoma
	// 		totalbuf += "\t";
	// 		totalbuf += "-"; // strand
	// 		totalbuf += "\t";
	// 		totalbuf += to_string(tmp_mismatch); // mismatch per target
	// 		totalbuf += "\t";
	// 		totalbuf += "0"; // num bulge, sempre uguale 0
	// 		totalbuf += "\t";
	// 		totalbuf += to_string(tmp_mismatch); // total per target
	// 		totalbuf += "\n";
	// 	}

	// 	// check delle posizioni con risultati negativi
	// 	for (int i = 0; i < resneg.size(); i++)
	// 	{
	// 		string target = reversetarget(genome.substr(resneg[i], guides[guideneg[i]].length()));
	// 		totalbuf += "X"; // bulge type (X stands for no bulge)
	// 		totalbuf += "\t";
	// 		totalbuf += guides[guideneg[i]];
	// 		totalbuf += "\t";
	// 		totalbuf += missmatching(target, resneg[i], guideneg[i], 1);
	// 		totalbuf += "\t";
	// 		totalbuf += chrnames;
	// 		totalbuf += "\t";
	// 		totalbuf += to_string(resneg[i]);
	// 		totalbuf += "\t";
	// 		totalbuf += to_string(resneg[i]); // Cluster Position
	// 		totalbuf += "\t";
	// 		totalbuf += "+";
	// 		totalbuf += "\t";
	// 		totalbuf += to_string(tmp_mismatch);
	// 		totalbuf += "\t";
	// 		totalbuf += "0";
	// 		totalbuf += "\t";
	// 		totalbuf += to_string(tmp_mismatch); // total per target
	// 		totalbuf += "\n";
	// 	}
	// }

	if (nowriting != 'r')
	{
		totalbuf.clear();
	}
	else
	{
		results << totalbuf;
		totalbuf.clear();
	}
}

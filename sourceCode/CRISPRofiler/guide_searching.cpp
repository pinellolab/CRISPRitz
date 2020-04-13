#include "include/crispritz.h"

using namespace std;

int tmp_mismatch;
string totalbuf;
extern ofstream results;
extern int guidelen, currentmissmatch, pamlimit, guidecount, totalguides, threads, guidelencorrected, inputmissmatch;
extern vector<int> missmatchthrestotal, pamindices, pamindicesreverse;
extern vector<string> guides;
extern string chrnames, genome;
extern char nowriting, noprofile;
extern bool pamdirection;

//profile
extern vector<vector<int>> guideprofiling;
extern vector<vector<vector<vector<int>>>> matrixprofiling;

//genome
extern vector<bitset<4>> genomebit;
extern vector<vector<bitset<4>>> guidesbit;
extern vector<vector<bitset<4>>> reverseguidesbit;

string reversetarget(string targetfound) //reverse del target per poterlo confrontare con la guida
{
#pragma omp parallel for
	for (int j = 0; j < guidelen; j++)
	{
		if (targetfound[j] == 'A')
		{
			targetfound[j] = 'T';
		}
		else if (targetfound[j] == 'T')
		{
			targetfound[j] = 'A';
		}
		else if (targetfound[j] == 'G')
		{
			targetfound[j] = 'C';
		}
		else if (targetfound[j] == 'C')
		{
			targetfound[j] = 'G';
		}
		else if (targetfound[j] == 'R')
		{
			targetfound[j] = 'Y';
		}
		else if (targetfound[j] == 'Y')
		{
			targetfound[j] = 'R';
		}
		else if (targetfound[j] == 'S')
		{
			targetfound[j] = 'S';
		}
		else if (targetfound[j] == 'W')
		{
			targetfound[j] = 'W';
		}
		else if (targetfound[j] == 'M')
		{
			targetfound[j] = 'K';
		}
		else if (targetfound[j] == 'K')
		{
			targetfound[j] = 'M';
		}
		else if (targetfound[j] == 'H')
		{
			targetfound[j] = 'D';
		}
		else if (targetfound[j] == 'D')
		{
			targetfound[j] = 'H';
		}
		else if (targetfound[j] == 'B')
		{
			targetfound[j] = 'V';
		}
		else if (targetfound[j] == 'V')
		{
			targetfound[j] = 'B';
		}
		else if (targetfound[j] == 'N')
		{
			targetfound[j] = 'N';
		}
	}

	reverse(targetfound.begin(), targetfound.end());
	return targetfound;
}

string missmatching(string targetfound, int index, int guidafound, int reverse) //conteggio dei mismatches e profilazione se richiesta dall'utente
{
	if (noprofile == 'p') //se richiesto profilo entro
	{
		//assegnazione valori ai contatori temporanei
		vector<int> tmp_mismatch_position;
		vector<vector<float>> tmp_bp_mismatch;
		tmp_mismatch = 0;
		tmp_mismatch_position.resize(guidelen);
		tmp_bp_mismatch.resize(guidelen);

		if (reverse == 0)
		{
			//inizio ricerca sulla guida
			for (int j = 0; j < guidelencorrected; j++)
			{
				tmp_bp_mismatch[j].resize(4, 0); //contatore per motiflogo

				if ((genomebit[index + j] & guidesbit[guidafound][j]) == 0)
				{
					int pos = 0;
					int pos_corrected = 0;
					if (!pamdirection)
					{
						pos = j;
						pos_corrected = pos;
					}
					else
					{
						pos = guidelen - 1 - j;
						pos_corrected = pos - pamlimit;
					}

					tmp_bp_mismatch[pos_corrected].resize(4, 0); //contatore per motiflogo

					//catena di elseif per verificare quale base genera MM
					if (targetfound[pos] == 'A')
						tmp_bp_mismatch[pos_corrected][0]++;
					else if (targetfound[pos] == 'T')
						tmp_bp_mismatch[pos_corrected][3]++;
					else if (targetfound[pos] == 'G')
						tmp_bp_mismatch[pos_corrected][2]++;
					else if (targetfound[pos] == 'C')
						tmp_bp_mismatch[pos_corrected][1]++;
					else if (targetfound[pos] == 'R')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + 0.5;
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + 0.5;
					}
					else if (targetfound[pos] == 'Y')
					{
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + 0.5;
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + 0.5;
					}
					else if (targetfound[pos] == 'S')
					{
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + 0.5;
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + 0.5;
					}
					else if (targetfound[pos] == 'W')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + 0.5;
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + 0.5;
					}
					else if (targetfound[pos] == 'K')
					{
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + 0.5;
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + 0.5;
					}
					else if (targetfound[pos] == 'M')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + 0.5;
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + 0.5;
					}
					else if (targetfound[pos] == 'B')
					{
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + (1 / 3);
					}
					else if (targetfound[pos] == 'D')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + (1 / 3);
					}
					else if (targetfound[pos] == 'H')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + (1 / 3);
					}
					else if (targetfound[pos] == 'V')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + (1 / 3);
					}

					targetfound[pos] = tolower(targetfound[pos]); //lettera minuscola
					tmp_mismatch++;										 //conto MM totali della guida
					if (!pamdirection)
					{
						guideprofiling[guidafound][pos]++; //segno un MM in una base della guida nel profilo
						tmp_mismatch_position[pos] = 5;	 //segno posizione MM
					}
					else
					{
						guideprofiling[guidafound][pos - pamlimit]++;
						tmp_mismatch_position[pos - pamlimit] = 5;
					}
				}
			}
		}
		else
		{
			for (int j = 0; j < guidelencorrected; j++)
			{
				tmp_bp_mismatch[j].resize(4, 0);

				if ((genomebit[index + j + pamlimit] & reverseguidesbit[guidafound][j]) == 0)
				{
					int pos;
					int pos_corrected = 0;
					if (!pamdirection)
					{
						pos = guidelencorrected - 1 - j;
						pos_corrected = pos;
					}
					else
					{
						pos = j + pamlimit;
						pos_corrected = pos - pamlimit;
					}
					tmp_bp_mismatch[pos_corrected].resize(4, 0);

					if (targetfound[pos] == 'A')
						tmp_bp_mismatch[pos_corrected][0]++;
					else if (targetfound[pos] == 'T')
						tmp_bp_mismatch[pos_corrected][3]++;
					else if (targetfound[pos] == 'G')
						tmp_bp_mismatch[pos_corrected][2]++;
					else if (targetfound[pos] == 'C')
						tmp_bp_mismatch[pos_corrected][1]++;
					else if (targetfound[pos] == 'R')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + 0.5;
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + 0.5;
					}
					else if (targetfound[pos] == 'Y')
					{
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + 0.5;
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + 0.5;
					}
					else if (targetfound[pos] == 'S')
					{
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + 0.5;
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + 0.5;
					}
					else if (targetfound[pos] == 'W')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + 0.5;
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + 0.5;
					}
					else if (targetfound[pos] == 'K')
					{
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + 0.5;
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + 0.5;
					}
					else if (targetfound[pos] == 'M')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + 0.5;
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + 0.5;
					}
					else if (targetfound[pos] == 'B')
					{
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + (1 / 3);
					}
					else if (targetfound[pos] == 'D')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + (1 / 3);
					}
					else if (targetfound[pos] == 'H')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][3] = tmp_bp_mismatch[pos_corrected][3] + (1 / 3);
					}
					else if (targetfound[pos] == 'V')
					{
						tmp_bp_mismatch[pos_corrected][0] = tmp_bp_mismatch[pos_corrected][0] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][1] = tmp_bp_mismatch[pos_corrected][1] + (1 / 3);
						tmp_bp_mismatch[pos_corrected][2] = tmp_bp_mismatch[pos_corrected][2] + (1 / 3);
					}

					targetfound[pos] = tolower(targetfound[pos]);
					tmp_mismatch++;
					if (!pamdirection)
					{
						guideprofiling[guidafound][pos]++;
						tmp_mismatch_position[pos] = 5;
					}
					else
					{
						guideprofiling[guidafound][pos - pamlimit]++;
						tmp_mismatch_position[pos - pamlimit] = 5;
					}
				}
			}
		}

		if (tmp_mismatch > 0)
		{
			//aggiorno il profili della guida con OFF-TARGET E MM TOTALI
			guideprofiling[guidafound][guidelencorrected + 1]++;
			guideprofiling[guidafound][guidelencorrected + 2] += tmp_mismatch;
		}
		else
		{
			//aggiorno profilo con ON-TARGET
			guideprofiling[guidafound][guidelencorrected]++;
		}

		//aggiorno profilo con nuovi MM nella colonna con threshold
		guideprofiling[guidafound][guidelen + tmp_mismatch]++;

		//aggiorno extendend_profile con dati della guida
		for (int gg = 0; gg < guidelencorrected; gg++)
		{
			if (tmp_mismatch_position[gg] == 5)
			{
				matrixprofiling[guidafound][tmp_mismatch][gg][4]++;
			}

			for (int dd = 0; dd < 4; dd++)
			{
				matrixprofiling[guidafound][tmp_mismatch][gg][dd] += tmp_bp_mismatch[gg][dd];
			}
		}
		return targetfound;
	}
	else if (nowriting == 'r' && noprofile != 'p') //entro in questo if se il profilo è disattivato, conto solo MM, metto minuscole le lettere che hanno generato MM e ritorno la guida
	{
		tmp_mismatch = 0;

		if (reverse == 0)
		{
			for (int j = 0; j < guidelencorrected; j++)
			{
				if ((genomebit[index + j] & guidesbit[guidafound][j]) == 0)
				{
					tmp_mismatch++;
					if (!pamdirection)
					{
						targetfound[j] = tolower(targetfound[j]);
					}
					else
					{
						targetfound[guidelen - 1 - j] = tolower(targetfound[guidelen - 1 - j]);
					}
				}
			}
		}
		else
		{
			for (int j = 0; j < guidelencorrected; j++)
			{
				if ((genomebit[index + j + pamlimit] & reverseguidesbit[guidafound][j]) == 0)
				{
					tmp_mismatch++;
					if (!pamdirection)
					{
						targetfound[guidelencorrected - 1 - j] = tolower(targetfound[guidelencorrected - 1 - j]);
					}
					else
					{
						targetfound[j + pamlimit] = tolower(targetfound[j + pamlimit]);
					}
				}
			}
		}
	}
	return targetfound;
}

//execute guide search
void guide_searching()
{
	vector<int> respos, resneg, respos_private, resneg_private;			  //inizializzazione variabili per store indici pam
	vector<int> guidepos, guideneg, guidepos_private, guideneg_private; //inizializzaione variabili per store indici guide

	int pampossize = pamindices.size();
	int pamnegsize = pamindicesreverse.size();

	//parallel region to find all targets on a selected chromosome
#pragma omp parallel num_threads(threads) private(respos_private, guidepos_private, resneg_private, guideneg_private)
	{
#pragma omp for simd schedule(static) nowait
		for (int i = 0; i < pampossize; i++) //ciclo sugli indici trovati precendetemente, PAM
		{
			for (int guidecount = 0; guidecount < totalguides; guidecount++) //ciclo le guide, su tutti gli indici
			{
				int currentmissmatch = 0; //missmatch di ogni guida

				for (int j = 0; j < guidelencorrected & (currentmissmatch <= inputmissmatch); j++) //confronto lettera-lettera, guida su genoma
				{
					if ((genomebit[j + pamindices[i]] & guidesbit[guidecount][j]) == 0)
					{
						currentmissmatch++;
					}
				}
				if (currentmissmatch <= inputmissmatch) //se rimango nella soglia, salvo dati, guida e pam che hanno generato target
				{
					respos_private.push_back(pamindices[i]);
					guidepos_private.push_back(guidecount);
				}
			}
		}

#pragma omp for simd schedule(static) nowait
		for (int i = 0; i < pamnegsize; i++)
		{
			for (int guidecount = 0; guidecount < totalguides; guidecount++)
			{
				int currentmissmatch = 0;

				for (int j = 0; j < guidelencorrected & (currentmissmatch <= inputmissmatch); j++)
				{
					if ((genomebit[j + pamlimit + pamindicesreverse[i]] & reverseguidesbit[guidecount][j]) == 0)
					{
						currentmissmatch++;
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
			//copy all threads values in a single value usable to write results
			respos.insert(respos.end(), respos_private.begin(), respos_private.end());
			guidepos.insert(guidepos.end(), guidepos_private.begin(), guidepos_private.end());

			resneg.insert(resneg.end(), resneg_private.begin(), resneg_private.end());
			guideneg.insert(guideneg.end(), guideneg_private.begin(), guideneg_private.end());
		}
	}

	if (!pamdirection) //direzione pam (se vero upstream altrimenti downstream)
	{
		//check delle posizioni con risultati positivi
		for (int i = 0; i < respos.size(); i++)
		{
			totalbuf += "X"; //bulge type (X stands for no bulge)
			totalbuf += "\t";
			totalbuf += guides[guidepos[i]]; //guida
			totalbuf += "\t";
			totalbuf += missmatching(genome.substr(respos[i], guidelen), respos[i], guidepos[i], 0); //calcolo della stringa e aggiornamento contatori profili
			totalbuf += "\t";
			totalbuf += chrnames; //nome chr
			totalbuf += "\t";
			totalbuf += to_string(respos[i]); //posizione del TARGET nel genoma
			totalbuf += "\t";
			totalbuf += to_string(respos[i]); //Cluster Position del TARGET nel genoma
			totalbuf += "\t";
			totalbuf += "+"; //strand
			totalbuf += "\t";
			totalbuf += to_string(tmp_mismatch); //mismatch per target
			totalbuf += "\t";
			totalbuf += "0"; //num bulge, sempre uguale 0
			totalbuf += "\t";
			totalbuf += to_string(tmp_mismatch); //total per target
			totalbuf += "\n";
		}

		//check delle posizioni con risultati negativi
		for (int i = 0; i < resneg.size(); i++)
		{
			totalbuf += "X"; //bulge type (X stands for no bulge)
			totalbuf += "\t";
			totalbuf += guides[guideneg[i]];
			totalbuf += "\t";
			totalbuf += missmatching(reversetarget(genome.substr(resneg[i], guidelen)), resneg[i], guideneg[i], 1);
			totalbuf += "\t";
			totalbuf += chrnames;
			totalbuf += "\t";
			totalbuf += to_string(resneg[i]);
			totalbuf += "\t";
			totalbuf += to_string(resneg[i]);	//Cluster Position
			totalbuf += "\t";
			totalbuf += "-";
			totalbuf += "\t";
			totalbuf += to_string(tmp_mismatch);
			totalbuf += "\t";
			totalbuf += "0";
			totalbuf += "\t";
			totalbuf += to_string(tmp_mismatch); //total per target
			totalbuf += "\n";
		}
	}
	else
	{
		for (int i = 0; i < respos.size(); i++)
		{
			totalbuf += "X"; //bulge type (X stands for no bulge)
			totalbuf += "\t";
			totalbuf += guides[guidepos[i]]; //guida
			totalbuf += "\t";
			totalbuf += missmatching(reversetarget(genome.substr(respos[i], guidelen)), respos[i], guidepos[i], 0); //calcolo della stringa e aggiornamento contatori profili
			totalbuf += "\t";
			totalbuf += chrnames; //nome chr
			totalbuf += "\t";
			totalbuf += to_string(respos[i]); //posizione del TARGET nel genoma
			totalbuf += "\t";
			totalbuf += to_string(respos[i]); //Cluster Position del TARGET nel genoma
			totalbuf += "\t";
			totalbuf += "-"; //strand
			totalbuf += "\t";
			totalbuf += to_string(tmp_mismatch); //mismatch per target
			totalbuf += "\t";
			totalbuf += "0"; //num bulge, sempre uguale 0
			totalbuf += "\t";
			totalbuf += to_string(tmp_mismatch); //total per target
			totalbuf += "\n";
		}

		//check delle posizioni con risultati negativi
		for (int i = 0; i < resneg.size(); i++)
		{
			totalbuf += "X"; //bulge type (X stands for no bulge)
			totalbuf += "\t";
			totalbuf += guides[guideneg[i]];
			totalbuf += "\t";
			totalbuf += missmatching(genome.substr(resneg[i], guidelen), resneg[i], guideneg[i], 1);
			totalbuf += "\t";
			totalbuf += chrnames;
			totalbuf += "\t";
			totalbuf += to_string(resneg[i]);
			totalbuf += "\t";
			totalbuf += to_string(resneg[i]);		//Cluster Position
			totalbuf += "\t";
			totalbuf += "+";
			totalbuf += "\t";
			totalbuf += to_string(tmp_mismatch);
			totalbuf += "\t";
			totalbuf += "0";
			totalbuf += "\t";
			totalbuf += to_string(tmp_mismatch); //total per target
			totalbuf += "\n";
		}
	}

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
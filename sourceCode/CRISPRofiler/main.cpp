#include "include/crispritz.h"

using namespace std;

ofstream results, profile, extentedprofile;
int threads, inputmissmatch;
double starttotal, endtotal;
string fastaname, pamname, guidename, resultname, profilename, extendedprofilename;
char nowriting, noprofile;
string writeprofile, writeextensiveprofile;
bool pamdirection, variant;

int main(int argc, char **argv)
{
   //start total execution time
   starttotal = omp_get_wtime();

   //assign argv variables to stream
   string resultwriting = "r";
   string profilewriting = "p";
   string profileplusresult = "t";
   nowriting = 'n';
   noprofile = 'n';
   fastaname = argv[1];
   pamname = argv[2];
   guidename = argv[3];
   inputmissmatch = atoi(argv[4]);
   resultname = argv[5];
   resultname += ".targets.txt";
   profilename = argv[5];
   profilename += ".profile.xls";
   extendedprofilename = argv[5];
   extendedprofilename += ".extended_profile.xls";
   pamdirection = 0;
   variant = atoi(argv[8]);

   //setting number threads used
   threads = atoi(argv[6]);
   if (threads > omp_get_max_threads() || threads == 0)
   {
      threads = omp_get_max_threads();
   }

   if (argc > 7 && (argv[7] == resultwriting)) //consenso a scrivere i result
   {
      nowriting = 'r';
      results.open(resultname);
      results << "#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\n";
   }
   else if (argc > 7 && (argv[7] == profilewriting)) //consenso a scrivere i profili (profile+extended_profile)
   {
      noprofile = 'p';
      profile.open(profilename);
      extentedprofile.open(extendedprofilename);
   }
   else if (argc > 7 && (argv[7] == profileplusresult)) //consenso a scrivere profili e result
   {
      nowriting = 'r';
      noprofile = 'p';
      results.open(resultname);
      profile.open(profilename);
      extentedprofile.open(extendedprofilename);
      results << "#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\n";
   }

   //reading pam
   reading_pam();

   //reading guides
   reading_guide();

   //profiler setting
   profilersetting();

   cout << "USED THREADS " << threads << endl;

   //reading chromosomes and execute analysis
   reading_chromosomes(argv); //inizio la ricerca sui cromosomi

   //profiling guides
   profiler(); //scrivo la profilazione

   //close the results file
   if (argc > 7 && (argv[7] == resultwriting))
   {
      results.close(); //chiudo il file result se era aperto
   }
   else if (argc > 7 && (argv[7] == profilewriting)) //chiudo i file di profili
   {
      profile << writeprofile;
      extentedprofile << writeextensiveprofile;

      profile.close();
      extentedprofile.close();
   }
   else if (argc > 7 && (argv[7] == profileplusresult)) //chiudo tutti i file aperti
   {
      results.close();

      profile << writeprofile;
      extentedprofile << writeextensiveprofile;
      profile.close();
      extentedprofile.close();
   }
   //end total execution time
   endtotal = omp_get_wtime();

   cout << "TOTAL TIME: " << endtotal - starttotal << endl;

   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "math.h"
#include "random.h"
#include <iomanip>
#include "funzioni.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


   rnd.SaveSeed();
 //____________________________________________________________//

  int M=10000;
  int N[4]={1,2,10,100};
  int N_int=100;

  int *v_std=new int[N_int];
  int *v_lor=new int[N_int];
  int *v_exp=new int[N_int];

  double S_std=0, S_lor=0, S_exp=0;
	int acc_std=0,acc_lor=0,acc_exp=0;
  double min_std=-1., max_std=1.;
  double min_lor=-10., max_lor=10.;
  double min_exp=0., max_exp=5.;

  double bin_std=(max_std-min_std)/double(N_int);
  double bin_lor=(max_lor-min_lor)/double(N_int);
  double bin_exp=(max_exp-min_exp)/double(N_int);


  ofstream output;
 
  string nomefile ="Histo_";
  string estensione=".txt";

  for(int k=0;k<4;k++){
	null(v_std,N_int);
	null(v_lor,N_int);
	null(v_exp,N_int);
	acc_std=0,acc_lor=0,acc_exp=0;
	for(int j=0;j<M;j++){
		S_std=0; S_lor=0; S_exp=0;		
		for(int l=0;l<N[k];l++){
			S_std+=rnd.Rannyu(-1.,1.);
			S_lor+=rnd.Cauchy_Lorentz(1,0);
			S_exp+=rnd.exponential(1);
			}

		S_std=S_std/double(N[k]);
		S_lor=S_lor/double(N[k]);
		S_exp=S_exp/double(N[k]);

		fill_Histo(v_std,N_int,min_std,bin_std,S_std,acc_std);
		fill_Histo(v_lor,N_int,min_lor,bin_lor,S_lor,acc_lor);
		fill_Histo(v_exp,N_int,min_exp,bin_exp,S_exp,acc_exp);
	}
	
	output.open(nomefile+to_string(N[k])+estensione);
	for(int s=0;s<N_int;s++){
		output<<setprecision(9)<<min_std+(s+1/2)*bin_std<<"\t"<<v_std[s]/(bin_std*acc_std)<<"\t"<<min_lor+(s+1/2)*bin_lor<<"\t"<<v_lor[s]/(bin_lor*acc_lor)<<"\t"<<min_exp+(s+1/2)*bin_exp<<"\t"<<v_exp[s]/(bin_exp*acc_exp)<<endl;
		}
		
	output.close();
  }

 delete[] v_std;
return 0;
}






/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

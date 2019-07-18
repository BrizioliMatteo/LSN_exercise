/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//Random numbers
#include "random.h"
#include "Chrom.h"
#include "Maps.h"
#include <vector>

int seed[4];
Random *rnd;

//parameters, observables
int N_p, N_c, rect, N_gen;

double r=2.; //esponente per la rigged roulet
double P_m=0.1;
double P_c=0.8;
double P_e=0.05;

std::vector<Chrom> pop;
//vector<Chrom> old_pop (N_p);
//vector<Chrom> new_pop (N_p);

Maps Mappa;

//pigreco

//functions
void Input(void); 
void Generation(void);
bool sort_fitness( Chrom &, Chrom &);
double Half_mean(void);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

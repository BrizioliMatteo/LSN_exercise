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

struct Var { double val;
				     int   rank; 
				  }; 

int seed[4];
Random *rnd;
Var Min_Cost;
//parameters, observables
int rect, N_c, nstep;
const int N_t=100;

double P_m=0.1;
double accepted,attempted;

Chrom new_male;
Chrom male;
Maps Mappa;


//geometric cooling
//pigreco

//functions
void Init_par(int);
void Input(int);
void Annealing(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

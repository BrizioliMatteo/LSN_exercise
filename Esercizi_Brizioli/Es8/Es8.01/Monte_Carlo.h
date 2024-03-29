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
int seed[4];
Random rnd;

//simulation
int nblk, nthrows, nstep;
const int m_props=2;
int n_props,ie,in;
const int nbins=100;
int conteggi=0;
// averages
double blk_norm,accepted,attempted;
double walker, blk_av;
double glob_av, glob_av2;
double histo[nbins];

double sigma;
double mu;
//configuration
double x;
double delta;
double bin_size,size;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void Measure(void);
void Equilibration(void);
double Error(double,double,int);
double Wave(double);
double Wave2(double);
double Ham(double);
void Save_histo(void);
void Fill_histo(void);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid_
#define __fluid_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iv, iw;
double vtail,ptail;
double walker[m_props];

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part];
std::string nomefile;
// thermodynamical state
int npart;
double beta,temp,vol,rho,box,rcut;

// simulation
int nstep;
double delta;
double accepted, attempted, acceptance_rate;
//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Instant_Value(int);
void Move(void);
void ConfFinal(void);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
void Equilibration(void);
#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

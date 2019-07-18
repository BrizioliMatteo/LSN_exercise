/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int M=10000;   //number of throw
const int steps=100; //number of Step
const int n_block=100;
const int L=int(M/n_block);
Random rnd;

double x[M]={0}, y[M]={0}, z[M]={0};
//pigreco
const double pi=3.1415927;
const double a=1.;//passo

//functions
void Input(void);
void Lattice_step(int);
void Continuum_step(int);
void Mean_position(int, const char *);
void null(double *, int);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

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
const int m_props=1000;
int iv,ik,it,ie,pr,n_props;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;
int Restart;

// averages
double glob_av[m_props], glob_av2[m_props], blk_av[m_props], walker[m_props];
double blk_norm;
//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;


//Valori numerici
//double sigma, e_kb, eps, m, delta;


// simulation
int nstep, iprint, seed, nblk;
double delta;
std::string nomefile;
//gdir
int nbins, igofr;
double bin_size;

const double pi=3.1415927;
//functions
void Input(void);
void Move(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double, double, int);
void Fill_Histo(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

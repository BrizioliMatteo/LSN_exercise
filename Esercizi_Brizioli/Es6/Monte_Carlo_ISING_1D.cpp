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
#include <ostream>
#include <cmath>
#include <string>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  if(Restart==1) SaveValue();
  ConfFinal(); //Write final configuration
  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadConf;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();
	
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> Restart;

	cout<<Restart<<endl;
  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

	
	 if(Restart==0){
   	ifstream input("seed.in");
   	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   	rnd.SetRandom(seed,p1,p2);
   	input.close();
	}else{
		ifstream input("seed.out");
   	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   	rnd.SetRandom(seed,p1,p2);
   	input.close();
	}


//Prepare arrays for measurements
  iu = 0; //Energy
  im = 1; //Heat capacity
  ix = 2; //Magnetic susceptibility
  ic = 3; //Magnetization
   
  n_props = 4; //Number of observables

  if(Restart==0){
    //initial configuration
    for (int i=0; i<nspin; ++i)
    {   
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }else if(Restart==1){
    ReadConf.open("Conf_final/config_T"+to_string(Dec(temp))+"_h"+to_string(Dec(h))+".final");
    for(int i=0;i<nspin; ++i)  ReadConf>>s[i];
    ReadConf.close();
  }


//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm, alpha;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);


    sm=s[o];
    energy_old=Boltzmann(sm,o);
    sm=-sm;//flip spin
    energy_new=Boltzmann(sm,o);


    if(metro==1) //Metropolis
    {
	p=exp(beta*(energy_old-energy_new));

	alpha=min(1.,p);
	if(alpha >= rnd.Rannyu()) {//Update
           s[o] =sm;
           accepted = accepted + 1.0;
         }
         attempted= attempted + 1.0;
    }else{  //Gibbs sampling 
    	p=1./(1.+exp(beta*(energy_old-energy_new)));   
		  if(p <= rnd.Rannyu())
	  	 s[o]=sm;
    }
  }

}




double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m +=s[i];
  }

 
  walker[iu] = u;
  walker[im] = m;
  walker[ix] = beta*m*m;
  walker[ic] = u*u;
  
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    cout << "Block number " << iblk << endl;
		if(metro==1) cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    ofstream InstConf;
    string nomefile;
    
    if(metro==1){
	 nomefile ="Metro_";
    }else{
	 nomefile ="Gibbs_";
    }


    InstConf.open("Instant_Conf/"+nomefile+"instant_T"+to_string(Dec(temp))+"_h"+to_string(Dec(h))+".0",ios::app);
   

    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    u_n=glob_av[iu]/(double)iblk;

     stima_x = blk_av[ix]/blk_norm/(double)nspin; //susceptibility
     glob_av[ix]  += stima_x;
     glob_av2[ix] += stima_x*stima_x;
     err_x=Error(glob_av[ix],glob_av2[ix],iblk);
     x_n=glob_av[ix]/(double)iblk;



     stima_c = beta*beta*(blk_av[ic]/blk_norm-(blk_av[iu]/blk_norm*blk_av[iu]/blk_norm))/(double)nspin; //Heat
     glob_av[ic]  += stima_c;
     glob_av2[ic] += stima_c*stima_c;
     err_c=Error(glob_av[ic],glob_av2[ic],iblk);
     c_n=glob_av[ic]/(double)iblk;
		
		if(h!=0.){
     stima_m = blk_av[im]/blk_norm/(double)nspin;
     glob_av[im]  += stima_m;
     glob_av2[im] += stima_m*stima_m;
     err_m=Error(glob_av[im],glob_av2[im],iblk);
     m_n=glob_av[im]/(double)iblk;
		 InstConf<<setprecision(9)<<iblk<<"\t"<<u_n<<"\t"<<c_n<<"\t"<<x_n<<"\t"<<m_n<<endl;
    }
		 else InstConf<<setprecision(9)<<iblk<<"\t"<<u_n<<"\t"<<c_n<<"\t"<<x_n<<endl;
    
		InstConf.close();
    cout << "----------------------------" << endl << endl;
}

void SaveValue(void){
   ofstream Ene, Heat, Mag, Chi;
   const int wd=9;
   string nomefile;
    if(metro==1){
	 nomefile ="Metro_";
    }else{
	 nomefile ="Gibbs_";
    }
    

    Ene.open("Result/"+nomefile+"ave.ene.dat",ios::app);
    Heat.open("Result/"+nomefile+"ave.heat.dat",ios::app);
    Mag.open("Result/"+nomefile+"ave.mag.dat",ios::app);
    Chi.open("Result/"+nomefile+"ave.sus.dat",ios::app);


    Ene<<setprecision(wd)<<temp<<"\t"<<h<<"\t"<<u_n<<"\t"<<err_u<<endl;
    Heat<<setprecision(wd)<<temp<<"\t"<<h<<"\t"<<c_n<<"\t"<<err_c<<endl;
    Chi<<setprecision(wd)<<temp<<"\t"<<h<<"\t"<<x_n<<"\t"<<err_x<<endl;
    if(h!=0.){
      Mag<<setprecision(wd)<<temp<<"\t"<<h<<"\t"<<m_n<<"\t"<<err_m<<endl;
    }

    Ene.close();
    Heat.close();
    Mag.close();
    Chi.close();
}


void ConfFinal(void)
{
  ofstream WriteConf;
  

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("Conf_final/config_T"+to_string(Dec(temp))+"_h"+to_string(Dec(h))+".final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}



int Dec(double var) 
{ 
    int value = (int)(var * 100); 
    return value;
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

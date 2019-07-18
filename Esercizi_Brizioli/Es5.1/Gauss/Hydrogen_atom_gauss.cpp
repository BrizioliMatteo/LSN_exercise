#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Hydrogen_atom_gauss.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  return 0;
}



void Input(void)
{
  ifstream ReadInput;

  cout << "Sampling a wave function of an Hydrogen atom" << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "The program uses Bohr units" << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
  ReadInput.open("input.dat"); //Read input
  ReadInput >> delta;
  ReadInput >> nthrows;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> wave;

  if(wave==0) nomefile="1s";
	else nomefile="2p";
  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

      x = rnd.Gauss(0.,delta);
      y = rnd.Gauss(0.,delta);
      z = rnd.Gauss(0.,delta);

    
    Equilibration();

  ofstream ConfXYZ, Distr;
  ConfXYZ.open("output.XYZ_"+nomefile+".dat");
	Distr.open("output.distr_"+nomefile+".dat");
	ConfXYZ.close();
	Distr.close();

  cout<<"Starting acceptance rate "<< accepted/attempted <<endl<<endl;



}


void Equilibration(void){

    double xold, yold, zold, xnew, ynew, znew, p, rnew, rold, acceptance_rate, alpha ;
    const int N_equ=500;
    accepted=0.;
    attempted=0.;
    for(int i=1;i<N_equ;i++){
     
     
      //Old
      xold = x;
      yold = y;
      zold = z;

      rold=sqrt(xold*xold+yold*yold+zold*zold);

      //New
      xnew = rnd.Gauss(xold,delta);
      ynew = rnd.Gauss(yold,delta);
      znew = rnd.Gauss(zold,delta);

      rnew=sqrt(xnew*xnew+ynew*ynew+znew*znew);

      //Metropolis test

      if(wave==0)	p= exp(2.*(rold-rnew));
      else        p= exp(rold-rnew)*znew/zold*znew/zold;

      alpha=min(1.,p);
    
      if(alpha >= rnd.Rannyu()){
    //Update
         x = xnew;
         y = ynew;
         z = znew;

         accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;
  }
   
  acceptance_rate=accepted/attempted;

   if (acceptance_rate<0.48){
	delta=delta+0.1;
	Equilibration();
    }else if(acceptance_rate>0.52){
	delta=delta-0.1;
	Equilibration();
   }

}










void Move(void){
  double xold, yold, zold, xnew, ynew, znew, p, rnew, rold, alpha;
   ofstream ConfXYZ;
   ConfXYZ.open("output.XYZ_"+nomefile+".dat",ios::app);
  //Old
    xold = x;
    yold = y;
    zold = z;

    rold=sqrt(xold*xold+yold*yold+zold*zold);

  //New
      xnew = rnd.Gauss(xold,delta);
      ynew = rnd.Gauss(yold,delta);
      znew = rnd.Gauss(zold,delta);

    rnew=sqrt(xnew*xnew+ynew*ynew+znew*znew);

  //Metropolis test
    if(wave==0)	p= exp(2.*(rold-rnew));
    else        p= exp(rold-rnew)*znew/zold*znew/zold;

    alpha=min(1.,p);
   
    if(alpha >= rnd.Rannyu()){
    //Update
       x = xnew;
       y = ynew;
       z = znew;

       ConfXYZ<<setprecision(9)<<x<<"\t"<<y<<"\t"<<z<<endl;

       accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;

    ConfXYZ.close();
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1){
           glob_av = 0;
           glob_av2 = 0;
   }

   blk_av = 0;
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{ 
   double r;
   r=sqrt(x*x+y*y+z*z);
   blk_av = blk_av + r;
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
   double r_mean,err;
   ofstream Distr;
    
    Distr.open("output.distr_"+nomefile+".dat",ios::app);
    
    r_mean = blk_av/blk_norm; //Potential energy
    glob_av += r_mean;
    glob_av2+= r_mean*r_mean;
    err=Error(glob_av,glob_av2,iblk);
    
    Distr<<setprecision(9)<<iblk<<"\t"<<glob_av/double(iblk)<<"\t"<<err<<endl;

    Distr.close();

}


double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}






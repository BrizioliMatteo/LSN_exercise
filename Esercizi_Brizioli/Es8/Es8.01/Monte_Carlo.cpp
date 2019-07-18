#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo.h"

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
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  Save_histo();
  return 0;
}



void Input(void)
{
  ifstream ReadInput;

  cout << "Monte Carlo simulation             " << endl << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
   


   //Read input informations
  ReadInput.open("input.dat");


  ReadInput >> mu;

  ReadInput >> sigma;

  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;
  


  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Delta = "<<delta<<endl;
  ReadInput.close();


  ofstream Inte, Out;//per evitare di sovrascrivere
  Out.open("Ave_int.dat");
  Inte.open("output.Integral.dat");
	Out.close();
	Inte.close();
  ie=0;
  in=1;
  n_props=2;

  x = rnd.Rannyu(-delta,delta);   

  Equilibration();

  cout<<"Starting acceptance rate "<< accepted/attempted <<endl<<endl;



}


void Equilibration(void){

    double xold, xnew, p, acceptance_rate, alpha ;
    const int N_equ=500;
    accepted=0.;
    attempted=0.;
    for(int i=1;i<N_equ;i++){
     
      //Old
      xold = x;

      //New
      xnew = rnd.Rannyu(-delta+xold,delta+xold);
      //Metropolis test
     p=Wave2(xnew)/Wave2(xold);
      alpha=min(1.,p);
    
      if(alpha >= rnd.Rannyu()){
    //Update
         x = xnew;
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
  double xold, xnew, p, alpha;
  //Old
    xold = x;

  //New
    xnew = rnd.Rannyu(-delta+xold,delta+xold);
    
    p=Wave2(xnew)/Wave2(xold);

  //Metropolis test
    alpha=min(1.,p);
   
    if(alpha >= rnd.Rannyu()){
    //Update
       x = xnew;
       accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1){
           glob_av = 0;
           glob_av2 = 0;
   }
   blk_av = 0.;
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{ 
   blk_av = blk_av + walker;

   blk_norm = blk_norm + 1.0;
}


void Measure(){
  walker=Ham(x)/Wave(x);
  Fill_histo();
}


double Ham(double t){
  double Vpot, Kin;
  Kin=-0.5*(-1./pow(sigma,2.)*Wave(t)+pow(t-mu,2.)/pow(sigma,4.)*exp(-0.5*pow((t-mu)/sigma,2.))+pow(t+mu,2.)/pow(sigma,4.)*exp(-0.5*pow((t+mu)/sigma,2.)));
  Vpot=(pow(t,4.)-5./2.*pow(t,2.))*Wave(t);
  return Kin+Vpot;
}


double Wave(double t){	
  double Eval;
	Eval=exp(-0.5*pow((t-mu)/sigma,2.))+exp(-0.5*pow((t+mu)/sigma,2.));
  return Eval;
}

double Wave2(double t){
  double Eval;
  Eval=Wave(t)*Wave(t);
  return Eval;
}


void Averages(int iblk) //Print results for current block
{
   double int_av,err;
   ofstream Inte, Out;
    
    Inte.open("output.Integral.dat",ios::app);
    
    int_av = blk_av/blk_norm; //Potential energy
    glob_av += int_av ;
    glob_av2+= int_av*int_av;
    err=Error(glob_av,glob_av2,iblk);
    
    Inte<<setprecision(9)<<iblk<<"\t"<<glob_av/double(iblk)<<"\t"<<err<<endl;
    if(iblk==nblk){
      Out.open("Ave_int.dat",ios::app);
      Out<<setprecision(9)<<"\t"<<glob_av/double(iblk)<<"\t"<<err<<"\t"<<mu<<"\t"<<sigma<<endl;
      Out.close();
    }
    Inte.close();

}


double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Fill_histo(void){
  size=8.;
  bin_size=size/double(nbins);
	int i=floor(x/bin_size)+floor(size/2./bin_size);
  if( i>=0 && i<100){
    conteggi++;
    histo[i]++;
	}
}

void Save_histo(void){
	ofstream Ist;
  std::string nomefile;
  nomefile=to_string((int)(mu*10.))+"_"+to_string((int)(sigma*10.));
  Ist.open("Distr/Prob_distr"+nomefile+".out");
  int wd=9;
  for(int i=0;i<nbins;i++)
    Ist<<setprecision(wd)<<i<<"\t"<<bin_size*(i+1./2.)-size/2.<<"\t"<<histo[i]/(double)conteggi/bin_size<<endl;
  Ist.close();
}





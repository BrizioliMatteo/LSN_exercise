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

int main()
{ 
  Input(); //Inizialization

  for(int i=0;i<N;i++){
	Integral();
	Integral_p();
	Block_mean(i);
	}


  cout<<int(2.3)<<" "<<int(0.6)<< " "<<int(1.2)<<endl;
  return 0;
}




void Input(void){

   //initialization random variable
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

}


void Integral(void) {
	double x;
	intg=0;
	for(int j=0;j<L;j++){
		x=rnd.Rannyu();
		intg+=pi/2.0*cos(pi/2.0*x);
	}
	intg=intg/L;
return;
}

void Integral_p(void) {
	double x;
	intg_p=0;
	for(int j=0;j<L;j++){
		x=rnd.pol();
		intg_p+=pi/4.0*cos(pi/2.0*x)/(1.-x);
	}
	intg_p=intg_p/L;
return;
}




void Block_mean(int i){

	ofstream output;
	output.open("integrale.txt",ios::app);

	sum+=intg;
	sum2+=intg*intg;

	sum_p+=intg_p;
	sum2_p+=intg_p*intg_p;

	double I_value=sum/double(i+1);
	double I_value_p=sum_p/double(i+1);
	double dev_std,dev_std_p;


	if(i==0){
		dev_std=0;
		dev_std_p=0;
	}else{
		dev_std=sqrt((sum2/double(i+1)-I_value*I_value)/double(i));
		dev_std_p=sqrt((sum2_p/double(i+1)-I_value_p*I_value_p)/double(i));
	}
	output<<setprecision(9)<<i+1<<"\t"<<I_value<<"\t"<<dev_std<<"\t"<<I_value_p<<"\t"<<dev_std_p<<endl;
	output.close();

return;
}






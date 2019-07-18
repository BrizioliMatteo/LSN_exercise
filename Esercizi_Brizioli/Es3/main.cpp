#include <iostream>
#include <fstream>
#include <string.h>
#include "math.h"
#include "random.h"
#include <iomanip>
#include "funzioni.h"

using namespace std;

int main(){ 
  Input(); //Inizialization Random

  ofstream output;
  output.open("Direct.txt");
 
  for(int i=0;i<n_blocks;i++){
	Direct_option();
	Block_mean(i,output);
	}
   output.close();
   sum_c=0., sum2_c=0., sum_p=0., sum2_p=0.;


  output.open("Discret.txt");

   for(int i=0;i<n_blocks;i++){
	Discrete_option();
	Block_mean(i,output);
	}
  output.close();
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
return;
}

void Direct_option(){
	P=0.;
	C=0.;
	double W;
	double S=0.;
	for(int j=0;j<L;j++){
		W=rnd.Gauss(0.,T);
		S=S_0*exp((r-0.5*sigma*sigma)*T+sigma*W);
		if (S-K>0)
			C+=exp(-r*T)*(S-K);
		else
			P+=exp(-r*T)*(K-S);
	}
	
	P=P/L;
	C=C/L;
}


void Discrete_option(){
	P=0.;
	C=0.;
	double Z;
	double S=0.;
	for(int j=0;j<L;j++){
		S=S_0;
		for(int i=0;i<100;i++){
			Z=rnd.Gauss(0.,1.);
			S=S*exp((r-0.5*sigma*sigma)*step_t+sigma*Z*sqrt(step_t));
		}
		if (S-K>0)
			C+=exp(-r*T)*(S-K);
		else
			P+=exp(-r*T)*(K-S);
	}
	
	P=P/L;
	C=C/L;
}


void Block_mean(int i,std::ofstream& output){
	

	sum_c+=C;
	sum2_c+=C*C;

	sum_p+=P;
	sum2_p+=P*P;

	double C_mean=sum_c/double(i+1);
	double P_mean=sum_p/double(i+1);
	double C_dev_std, P_dev_std;

	if(i==0){
		C_dev_std=0.;
		P_dev_std=0.;
	}else{
		C_dev_std=sqrt((sum2_c/double(i+1)-C_mean*C_mean)/double(i));
		P_dev_std=sqrt((sum2_p/double(i+1)-P_mean*P_mean)/double(i));
	}
output<<setprecision(9)<<i+1<<"\t"<<C_mean<<"\t"<<C_dev_std<<"\t"<<P_mean<<"\t"<<P_dev_std<<endl;


return;
}



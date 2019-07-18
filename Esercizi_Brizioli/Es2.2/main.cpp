#include <iostream>
#include <fstream>
#include <string>
#include "math.h"
#include "random.h"
#include <iomanip>
#include "funzioni.h"

using namespace std;

int main(){ 
  Input(); //Inizialization Random
	const char nomefile[]="Lattice_random_walk.txt";
	const char nomefile2[]="Uniform_random_walk.txt";
  	for(int j=0;j<steps;j++){
		for(int i=0;i<M;i++){
	  		Lattice_step(i);
		}
		Mean_position(j,nomefile);
	}	
	
	null(x,M);
	null(y,M);
	null(z,M);
	
  	for(int j=0;j<steps;j++){
		for(int i=0;i<M;i++){
	  		Continuum_step(i);
		}
		Mean_position(j,nomefile2);
	}	

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







void Lattice_step(int i){
	int r= int(rnd.Rannyu(0.,3.));	
	int s=int(rnd.Rannyu(0.,2.));
	if (s==0)
		s=-1;
	
	if (r==0)
		x[i]=x[i]+s;
	else if(r==1)
		y[i]=y[i]+s;
	else
		z[i]=z[i]+s;

return;
}

void Continuum_step(int i){
	double theta = rnd.Theta();	
	double phy   = rnd.Rannyu(0.,2.*pi);

	x[i]+=sin(theta)*cos(phy)*a;
	y[i]+=sin(theta)*sin(phy)*a;
	z[i]+=cos(theta)*a;

return;

}

void Mean_position(int i,const char *nomefile){
	int k=0;
	double r_mean_block=0;
	double sum=0, sum2=0, mean_value, dev_std;
	for(int j=1;j<n_block;j++){
		for(int t=1;t<L;t++){
			k=t+L*j;
			r_mean_block+=double(x[k]*x[k]+y[k]*y[k]+z[k]*z[k]);
			
		}
		r_mean_block=r_mean_block/double(L);
		sum+=r_mean_block;
		sum2+=r_mean_block*r_mean_block;
	}

	mean_value=sum/(double (n_block));
	dev_std=sqrt((sum2/double(n_block)-mean_value*mean_value)/double(n_block-1));

	ofstream output;
	output.open(nomefile,ios::app);
	output<<setprecision(9)<<i+1<<"\t"<<sqrt(mean_value)<<"\t"<<sqrt(dev_std)<<endl;
	output.close();

return;
}



void null(double * v, int n){
	for(int j=0;j<n;j++)
		v[j]=0;
return;
}


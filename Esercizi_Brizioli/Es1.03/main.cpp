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
 
int main (int argc, char *argv[]){

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
 //____________________________________________________________//

 


	for(int i=0;i<N;i++){
		Media_blocco(); //media sul singolo blocco
		block(i);
	}

return 0;
}










void Media_blocco(){
	double y, y_proj,ang, dist;
	int N_bar;
	int N_hit=0;

	for(int i=0;i<L;i++){
 		y=rnd.Rannyu();    //genero la y e e l'angolo
 		ang=rnd.Theta();
		
		N_bar=floor(y/d);
 		y_proj=l*sin(ang);
		 if (y_proj<=0){	
		    	dist=abs(y-N_bar*d);
			//cout<<dist<<"\t"<<abs(y_proj)<<endl;	
			if (abs(y_proj)>dist)
				N_hit++;
		 }else if (y_proj>0){
		    dist=abs(y-(N_bar+1)*d);;	
			if (abs(y_proj)>dist)
				N_hit++;
		 }
	}

        media=2.*l*double(L)/(d*double(N_hit));

return;
}


void block(int i){

	ofstream output;
	output.open("PI.txt",ios::app);

	sum+=media;
	sum2+=media*media;

	double PI=sum/double(i+1);
	double dev_std;

	if(i==0){
		dev_std=0;
	}else{
		dev_std=sqrt((sum2/double(i+1)-PI*PI)/i);
	}

	output<<setprecision(9)<<i<<"\t"<<PI<<"\t"<<dev_std<<endl;
	output.close();

return;
	
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

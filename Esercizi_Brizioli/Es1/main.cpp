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


using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
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
	
	ofstream output,output2;
	
	output.open("average.out");
	output2.open("sigma_squared.out");

	int M=200000;
	int N=100;
	int L=M/N;
		
	double r=0;	

	double sum=0;
	double sum_sigma=0;

	double sum1=0;
	double sum2=0;
	double sum3=0;
	double sum4=0;
	
	double mean=0;
	double dev_std=0;
	double sigma=0;
	double dev_std_sigma=0;		

	
	for(int i=0;i<N;i++){
		sum=0;
		sum_sigma=0;
		for(int j=0;j<L;j++){
			r=rnd.Rannyu();	
			sum+=r;
			sum_sigma+=(r-0.5)*(r-0.5);
		}

		sum1+=sum/L;       //sommo le medie per ogni blocco
		sum2+=sum/L*sum/L;

		sum3+=sum_sigma/L;
		sum4+=sum_sigma/L*sum_sigma/L;	    //sommo le medie del quadrato per ogni blocco

		mean=sum1/(i+1);
		sigma=sum3/(i+1);
		if(i==0){
			dev_std=0;
			dev_std_sigma=0;
		}else{
			dev_std=sqrt((sum2/(i+1)-mean*mean)/i);
			dev_std_sigma=sqrt((sum4/(i+1)-sigma*sigma)/i);
		}
		output<<setprecision(9)<<i<<"	"<<mean<<"	"<<dev_std<<endl;
		output2<<setprecision(9)<<i<<"	"<<sigma<<"	"<<dev_std_sigma<<endl;
		
	}

	output.close();
	output2.close();



//____________________________________________________________//


	/*double sigma=0;
	double dev_std_sigma=0;	
	double r=0;
	sum1=0;
	sum2=0;
	

	for(int i=0;i<N;i++){
		sum=0;
		for(int j=0;j<L;j++){
			r=rnd.Rannyu();
			sum+=(r-0.5)*(r-0.5);
			}
		
		sum1+=sum/L;             //sommo le medie per ogni blocco
		sum2+=sum/L*sum/L;	    //sommo le medie del quadrato per ogni blocco
	
		sigma=sum1/(i+1);
		if(i==0)
			dev_std=0;
		else 
			dev_std_sigma=sqrt((sum2/(i+1)-sigma*sigma)/i);
		
		cout<<sigma<<" +- "<<dev_std_sigma<<endl;
		
	}*/




	



	
//____________________________________________________________//

	output.open("chi_squared.out");	
	
	int Int=100;
	int n=10000;
  	int count[Int]={0};
	double chi_squared=0;
	int i=0;

	for(int j=0;j<100;j++){
		for(int l=0;l<n;l++){
			i=int(rnd.Rannyu()*Int);
			count[i]++;	
		}
		chi_squared=0;
		for(int m=0;m<Int;m++){
		chi_squared+=((double(count[m])-double(n)/double(Int))*(double(count[m])-double(n)/double(Int)))/(double(n)/double(Int));
		count[m]=0;
		}
		output<<setprecision(9)<<j<<"	"<<chi_squared<<endl;
	}




 	output.close();










return 0;
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

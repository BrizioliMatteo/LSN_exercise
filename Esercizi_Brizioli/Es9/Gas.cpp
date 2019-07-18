#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include<algorithm>
#include "Gas.h" 


using namespace std;

int main(){
  Input(); //Inizialization
	ofstream Best, Half;
	Best.open("Best.dat");
	Half.open("Half.dat");
	for(int i=0;i<N_gen;i++){
  	Generation();
			cout<<"Generazione :"<<"\t"<<i<<endl;
			cout<<"Lunghezza percorso :"<<"\t"<<pop[0].Cost_fn_L1(Mappa)<<endl;
			cout<<"------------------------------------"<<endl;
			Best<<setprecision(9)<<i<<"\t"<<pop[0].Cost_fn_L1(Mappa)<<endl;
			Half<<setprecision(9)<<i<<"\t"<<Half_mean()<<endl;
	}
	pop[0].Save_conf(Mappa);
	Half.close();
	Best.close();
  delete rnd;
  return 0;
}




void Input(void){
  
  ifstream ReadInput;

//Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  rnd =new Random;
  
  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3]; 
  rnd->SetRandom(seed,p1,p2);
   
  input.close();   
  //Read input informations
  ReadInput.open("input.dat");


  ReadInput >> N_p;

  ReadInput >> N_c;

  ReadInput >> N_gen;

  ReadInput >> rect;


  cout << "Salesman problem" << endl;
  cout << "Numero di città da visitare= " <<N_c << endl;
  cout << "Numero di cromosomi nella popolazione= " << N_p << endl;
  cout << "Numero di generazioni = " << N_gen << endl;
  cout << "Distribuzione città = "<<rect<<endl;

  ReadInput.close();
 
	
  Mappa.Set_par(N_c,rect,rnd);
	
  for(int i=0;i<N_p;i++)		pop.push_back(Chrom(N_c,rnd));
  sort(pop.begin(), pop.end(), sort_fitness);

	}


void Generation(void){
  int m,f,cross;
	vector<Chrom> pop_new;
  for(int i=0;i<N_p/2;i++){//selection
		m= int(pow(rnd->Rannyu(0.,1.),r)*N_p);
    f= int(pow(rnd->Rannyu(0.,1.),r)*N_p);
		while(m==f) f= int(pow(rnd->Rannyu(0.,1.),r)*N_p);
    Chrom parent_m=pop[m];
    Chrom parent_f=pop[f];
		
		if(P_c>rnd->Rannyu()){//crossover
			cross=(int)rnd->Rannyu(0.,N_c-1);
			pop_new.push_back(Chrom(parent_m,parent_f,cross,rnd));
			pop_new.push_back(Chrom(parent_f,parent_m,cross,rnd));	
		}else{
			pop_new.push_back(parent_m);
			pop_new.push_back(parent_f);
		}
	}

	for(int i=0;i<N_p;i++){ 
		pop_new[i].Pair_permutation();
		pop_new[i].Shift();
		pop_new[i].Shift_m();
		pop_new[i].Inversion();
	}
	sort(pop_new.begin(), pop_new.end(), sort_fitness);
	
	if(P_e>rnd->Rannyu()){//elite
		pop_new.back()=pop.front();
		sort(pop_new.begin(), pop_new.end(), sort_fitness);
	}

	copy(pop_new.begin(),pop_new.end(),pop.begin());
}


double Half_mean(void){
	double mean_val=0;
	for(int i=0;i<N_p/2;i++){
		mean_val+=pop[i].Cost_fn_L1(Mappa);
	}
	return mean_val/double(N_p/2);
}



bool sort_fitness( Chrom &former, Chrom &latter){
	return former.Cost_fn_L1(Mappa)<latter.Cost_fn_L1(Mappa);
}






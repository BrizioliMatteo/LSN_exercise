#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include<algorithm>
#include "Gas.h" 
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);

	int my_rank,size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	Var local_loss;
	
	ofstream Best;


	Init_par(my_rank);
	Input(my_rank);
	Best.open("Data/Best_"+to_string(my_rank)+"_"+to_string(rect)+".dat");
	double beta=1.;
	for(int i=0;i<N_t;i++){
  	Annealing(beta);
		Best<<setprecision(9)<<i<<"\t"<<male.Cost_fn_L1(Mappa)<<endl;
		beta=beta*1.1;	
	}
	Best.close();
	
	local_loss.val=male.Cost_fn_L1(Mappa);
	local_loss.rank=my_rank;
	cout<<"Rank :"<<"\t"<<local_loss.rank<<endl;	
	cout<<"Lunghezza percorso :"<<"\t"<<local_loss.val<<endl;
	cout<<"------------------------------------"<<endl;
	

	MPI_Allreduce(&local_loss,&Min_Cost,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
	if(Min_Cost.rank==my_rank){
		cout<<"Processo che minimizza la distanza"<<endl;
		cout<<Min_Cost.val<<"\t"<<Min_Cost.rank<<endl;
		male.Save_conf(Mappa);
	}

	delete rnd;
	MPI_Finalize();
  return 0;
}

//con lo stesaso tempo di esecuzione ho più statistica in base al numero dei processi 


void Input(int my_rank){
	ifstream Read;
	Read.open("Primes");//Read input informations
	if(!Read.is_open()){
		cerr <<"Error: file Primes does not exist" <<endl;
		exit(-1);
	}	

	int p1,p2,count=0;
	while(count<=my_rank){
		Read >> p1 >> p2;
		count++;
	}
	Read.close();

	rnd->SetRandom(seed,p1,p2);
	male=Chrom(N_c,rnd);
}

void Init_par(int my_rank){
	int p1,p2;
	ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
	
  rnd =new Random;

	if(my_rank==0){
		ifstream input("seed.in");
 		input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  	input.close();


		ifstream ReadInput;
  	ReadInput.open("input.dat");

 		ReadInput >> N_c;
		ReadInput >> nstep;
  	ReadInput >> rect;
		
		ReadInput.close();
	}
	MPI_Bcast(&seed,4,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Bcast(&N_c,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nstep,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&rect,1,MPI_INT,0,MPI_COMM_WORLD);
		

  cout << "Salesman problem" << endl;
  cout << "Numero di città da visitare= " <<N_c << endl;
  cout << "Numero di step per temperatura = " << nstep << endl;
  cout << "Distribuzione città = "<<rect<<endl;

	rnd->SetRandom(seed,p1,p2);
  Mappa.Set_par(N_c,rect,rnd);
}



void Annealing(double beta){
	double p,alpha;
	accepted=0.;	
	attempted=0.;

	for(int i=0;i<nstep;i++){
		new_male=male;
	
		new_male.Pair_permutation();
		new_male.Shift();
		new_male.Shift_m();
		new_male.Inversion();

  	p=exp(beta*(male.Cost_fn_L1(Mappa)-new_male.Cost_fn_L1(Mappa)));

		alpha=min(1.,p);
		if(alpha >= rnd->Rannyu()) {	
   	 male=new_male;
   	 accepted = accepted + 1.0;
   	}
  	attempted= attempted + 1.0;
	}
}







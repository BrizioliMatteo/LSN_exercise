#include<iostream>
#include<algorithm>
#include"math.h"
#include"Chrom.h"
#include <fstream>
#include <ostream>
#include <iomanip>
using namespace std;


Chrom::Chrom(int N_all, Random *Srnd){
  m_N=N_all;
  Trnd=Srnd;
  for(int j=0;j<N_all;j++) vect.push_back(j);

  random_shuffle(vect.begin(),vect.end()); 
}
      

Chrom::~Chrom(){}

double Chrom::Cost_fn_L1(const Maps Cities){	
  double L1=Cities.Dist(vect.at(0),vect.at(m_N-1));
  for (unsigned int i=0; i<vect.size()-1; i++) L1+=Cities.Dist(vect.at(i),vect.at(i+1));
	return L1;
}



void Chrom::Pair_permutation(void){
  if(P_m>Trnd->Rannyu()){
 	 int ind_1=int(Trnd->Rannyu(0.,m_N));
 	 int ind_2=int(Trnd->Rannyu(0.,m_N));
 	 while(ind_1==ind_2)	ind_2=int(Trnd->Rannyu(0.,m_N));
 	 swap(vect[ind_1],vect[ind_2]);
  }
}


void Chrom::Shift(void){
  if(P_m>Trnd->Rannyu()){
 	 int ind=int(Trnd->Rannyu(1.,m_N-1));
   rotate(vect.begin(),vect.end()-ind,vect.end());
	}
}

void Chrom::Print_conf(void){
  for(int i=0;i<m_N;i++)
	  cout<<vect[i]<<" ";

  cout<<"\n";

}

void Chrom::Shift_m(void){
   
  if(P_m>Trnd->Rannyu()){//scelgo due intervalli di m contigui
		int m=int(Trnd->Rannyu(0.,m_N/2.));
		int n=int(Trnd->Rannyu(m_N/2,m_N-m));
    int ind=int(Trnd->Rannyu(0.,n-m));
    vector<int> vect_copy (m_N);
	
    copy(vect.begin()+ind,vect.begin()+ind+m,vect_copy.begin()+n);
    copy(vect.begin()+n,vect.begin()+n+m,vect_copy.begin()+ind);
    swap_ranges(vect_copy.begin()+ind,vect_copy.begin()+ind+m,vect.begin()+ind);
    swap_ranges(vect_copy.begin()+n,vect_copy.begin()+n+m,vect.begin()+n);
	}
}


void Chrom::Inversion(void){

  if(P_m>Trnd->Rannyu()){
   int m=int(Trnd->Rannyu(0.,m_N));
   int ind=int(Trnd->Rannyu(0.,m_N-m));
	 reverse(vect.begin()+ind,vect.end()-m);
	}
}

int Chrom::Get_Size(void)const{
	return vect.size();
}

int Chrom::Get_Element(int i)const{
	return vect.at(i);
}


Chrom::Chrom(const Chrom& parent1, const Chrom& parent2, int cross, Random *Srnd){
	m_N=parent1.Get_Size();
	Trnd=Srnd;
	for(int i=0;i<cross;i++) vect.push_back(parent1.Get_Element(i)); 
	int j=0;
	 while((int)vect.size()<m_N){
  			if(std::find(vect.begin(),vect.end(),parent2.Get_Element(j)) == vect.end())  vect.push_back(parent2.Get_Element(j));
				j++;
	}
}
 
void Chrom::Save_conf(const Maps Cities){
	ofstream Conf_cit("Maps");
  for(int i=0;i<m_N;i++)	Conf_cit<<setprecision(9)<<i<<"\t"<<Cities.GetPos(vect.at(i)).x<<"\t"<<Cities.GetPos(vect.at(i)).y<<endl;
	Conf_cit<<setprecision(9)<<m_N<<"\t"<<Cities.GetPos(vect.at(0)).x<<"\t"<<Cities.GetPos(vect.at(0)).y<<endl;	Conf_cit.close();

}

/*void Chrom::Check(void){//oppure copio il vettore faccio il sort e poi confronto a due a due
	for(int i=0;i<int(vect.size())-1;i++){	
		if (find(vect.begin()+i+1, vect.end(), vect[i]) != vect.end()){
				cout << "individual doesn't fulfil the bonds"<<endl;
				break;
		}
	}
}
  */





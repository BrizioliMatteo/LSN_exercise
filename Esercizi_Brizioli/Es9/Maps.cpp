#include"Maps.h"
#include<iostream>
#include<stdexcept>
#include<stdlib.h>
#include<sstream>
#include"math.h"

using namespace std;

Maps::Maps(){}





Maps::Maps(int Npts, int disp, Random *Trnd){
  m_N=Npts;
	Cities.resize(m_N);
 double theta;
	for(int j=0;j<m_N; j++){
    if(disp==0){//i punti sono generati su un quadrato lato 1
      Cities[j].x=Trnd->Rannyu(-1.,1.);
      Cities[j].y=Trnd->Rannyu(-1.,1.);
    }else{
      theta=Trnd->Rannyu(0.,2.*pi);
      Cities[j].x=cos(theta);
      Cities[j].y=sin(theta);
		}
	}
}



void Maps::Set_par(int Npts, int disp, Random *Trnd){
  m_N=Npts;
  
	Cities.resize(m_N);

 double theta;
	for(int j=0;j<m_N; j++){
    if(disp==0){//i punti sono generati su un quadrato lato 1
      Cities[j].x=Trnd->Rannyu(-1.,1.);
      Cities[j].y=Trnd->Rannyu(-1.,1.);
    }else{
      theta=Trnd->Rannyu(0.,2.*pi);
      Cities[j].x=cos(theta);
      Cities[j].y=sin(theta);
		}
	}
}
      

Pos Maps::GetPos(int i) const{
	if (i<m_N) return Cities[i];
	else{
	 cerr<<"l'intero assegnato e' maggiore della lunghezza della struct"<<endl;
	exit(-1);
	}
}


double Maps::Dist(int i, int l) const{
	double ds =pow(Cities[i].x-Cities[l].x,2.)+pow(Cities[i].y-Cities[l].y,2.);
  return sqrt(ds);
}



Maps::~Maps(){}


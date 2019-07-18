#include "funzioni.h"

void fill_Histo(int * v, int Int, double min, double bin, double S, int &accepted){
	int i=floor(S/bin)-floor(min/bin);
	if(i<Int && i>=0){
		accepted++;
		v[i]++;
	}
}

void null(int * v, int n){
	for(int j=0;j<n;j++)
		v[j]=0;
}

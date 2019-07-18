/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>         // rint, pow
#include <iomanip>
#include <string>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();       //Inizialization
  int nconf = 1;
   for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep/nblk; ++istep){
      Move();           //Move particles with Verlet algorithm
      if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
      if(istep%10 == 0){
        Measure();//Properties measurement
	Accumulate();
//      ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk);    //Print results for current block
  }
  ConfFinal();         //Write final configuration to restart
  return 0;
}












void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double fx[m_part],fy[m_part],fz[m_part];
  double xnew, ynew, znew;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> nblk;
  ReadInput >> Restart;  

  if(rho==1.1) nomefile="solid";
  else if(rho==0.8) nomefile="liquid";
  else if(rho==0.05) nomefile="gas";

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  pr = 4; //Pression
  n_props = 5; //Number of observables //Number of observables
//Read initial configuration
  igofr=5;
  nbins=100;
  n_props=n_props+nbins;
  bin_size = (box/2.0)/(double)nbins;

  if (Restart==0){
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
	    ReadConf >> x[i] >> y[i] >> z[i];
	    x[i] = x[i] * box;
	    y[i] = y[i] * box;
	    z[i] = z[i] * box;
	  }
	  ReadConf.close();

	//Prepare initial velocities
	  cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	  double sumv[3] = {0.0, 0.0, 0.0};
	  for (int i=0; i<npart; ++i){
	   vx[i] = rand() - 0.5;
	   vy[i] = rand() - 0.5;
	   vz[i] = rand() - 0.5;

	   sumv[0] += vx[i];
	   sumv[1] += vy[i];
	   sumv[2] += vz[i];
	  }
	  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	  double sumv2 = 0.0, fs;
	  for (int i=0; i<npart; ++i){
	    vx[i] = vx[i] - sumv[0];
	    vy[i] = vy[i] - sumv[1];
	    vz[i] = vz[i] - sumv[2];

	    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	  }
	  sumv2 /= (double)npart;

	  fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	  for (int i=0; i<npart; ++i){
	    vx[i] *= fs;
	    vy[i] *= fs;
	    vz[i] *= fs;

	    xold[i] = Pbc(x[i] - vx[i] * delta);
	    yold[i] = Pbc(y[i] - vy[i] * delta);
	    zold[i] = Pbc(z[i] - vz[i] * delta);
	  }
  }else{
	   //read r⃗ (t)and r⃗ (t−dt)r→(t−dt) from the corresponding input files
	   cout << "Read actual configuration from file config.final " << endl << endl;
	   ReadConf.open("Conf_final/config_"+nomefile+".final");
	   for (int i=0; i<npart; ++i){
	     ReadConf >> x[i] >> y[i] >> z[i];
	     x[i] = x[i] * box;
	     y[i] = y[i] * box;
	     z[i] = z[i] * box;
	   }
	   ReadConf.close(); 

	   cout << "Read old configuration from file old.final " << endl << endl;
	   ReadConf.open("Conf_final/old_"+nomefile+".final");
	   for (int i=0; i<npart; ++i){
	     ReadConf >> xold[i] >> yold[i] >> zold[i];
	     xold[i] = xold[i] * box;
	     yold[i] = yold[i] * box;
	     zold[i] = zold[i] * box;
	   }
     ReadConf.close(); 


           
	   //compute r⃗ (t+dt) with one step of the Verlet algorithm; with r⃗ (t+dt)and r⃗ (t) compute v⃗ (t+dt/2); finally obtain the actual temperature T(t+dt/2)
	   for(int i=0; i<npart; ++i){ //Force acting on particle i
	     fx[i] = Force(i,0);
	     fy[i] = Force(i,1);
	     fz[i] = Force(i,2);

	   }  

	   double sumv2 = 0.0, fs;
	   
	   for(int i=0; i<npart; ++i){ //Verlet integration scheme
             
	     xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
	     ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
	     znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

	     vx[i] = Pbc(xnew - x[i])/( delta);
	     vy[i] = Pbc(ynew - y[i])/( delta);
	     vz[i] = Pbc(znew - z[i])/( delta);

	     sumv2+= vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];

	     x[i] = xnew;
	     y[i] = ynew;
	     z[i] = znew;

	   }

	   //by comparing T(t+dt/2) with the desired/target temperature T⋆ extract a scaling factor for the velocities and rescale them: v⃗ (t)
	   sumv2 /= (double)npart;
           fs = sqrt(3. * temp / sumv2);   // fs = velocity scale factor 

     cout<<"Target Temperature "<<temp<<endl;
     cout<<"Actual Temperature "<<sumv2/3.<<endl;
 	   cout<<"Scale factor"<<fs<<endl;

	   for (int i=0; i<npart; ++i){
	     vx[i] *= fs;
	     vy[i] *= fs;
	     vz[i] *= fs;
	     xold[i] =Pbc(x[i] - vx[i] * delta);
	     yold[i] =Pbc(y[i] - vy[i] * delta);
	     zold[i] =Pbc(z[i] - vz[i] * delta);
	   }
          
   }
    return;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0.;
           glob_av2[i] = 0.;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0.;
   }
   
   blk_norm = 0.;
   return;
}






void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}



void Averages(int iblk) //Print results for current block
{
   double errore[m_props];
   double mean, r, delta_v,gdir;
   ofstream Epot, Ekin, Etot, Temp, Pres, Gofr, Gave;
   int wd=9;
   cout << "Block number " << iblk << endl;

   Epot.open("Output/ave_epot_"+nomefile+".out",ios::app);
   Ekin.open("Output/ave_ekin_"+nomefile+".out",ios::app);
   Etot.open("Output/ave_etot_"+nomefile+".out",ios::app);
   Temp.open("Output/ave_temp_"+nomefile+".out",ios::app);
   Pres.open("Output/ave_pres_"+nomefile+".out",ios::app);
   Gofr.open("Output/ave_gfor_"+nomefile+".out",ios::app);
   Gave.open("Gofr/Gave_"+nomefile+".out");

   for(int i=0;i<igofr;i++){
     mean=blk_av[i]/blk_norm;
     glob_av[i]+=mean;
     glob_av2[i]+=mean*mean;
     errore[i]=Error(glob_av[i],glob_av2[i],iblk);      
   }
 
   for (int k=igofr; k<igofr+nbins; ++k){
      r=bin_size*double(k-igofr);
      delta_v=rho*m_part*4./3.*pi*(pow(r+bin_size,3.)-pow(r,3.));
      gdir=blk_av[k]/blk_norm/delta_v;
      glob_av[k] += gdir;
      glob_av2[k] +=gdir*gdir;
      errore[k]=Error(glob_av[k],glob_av2[k],iblk);

      Gofr<<setprecision(wd)<<iblk<<"\t"<<r<<"\t"<<gdir<<"\t"<<glob_av[k]/double(iblk)<<"\t"<<errore[k]<<endl;
      if(iblk==nblk){
          Gave<<setprecision(wd)<<r<<"\t"<<glob_av[k]/double(iblk)<<"\t"<<errore[k]<<endl;
      }
    }

   

    
   Epot<<setprecision(wd)<<iblk<<"\t"<<glob_av[iv]/double(iblk)<<"\t"<<errore[iv]<<endl;
   Ekin<<setprecision(wd)<<iblk<<"\t"<<glob_av[ik]/double(iblk)<<"\t"<<errore[ik]<<endl;
   Etot<<setprecision(wd)<<iblk<<"\t"<<glob_av[ie]/double(iblk)<<"\t"<<errore[ie]<<endl;
   Temp<<setprecision(wd)<<iblk<<"\t"<<glob_av[it]/double(iblk)<<"\t"<<errore[it]<<endl;
   Pres<<setprecision(wd)<<iblk<<"\t"<<glob_av[pr]/double(iblk)<<"\t"<<errore[pr]<<endl;

   cout << "----------------------------" << endl << endl;

   Epot.close();
   Ekin.close();
   Temp.close();
   Etot.close();
   Pres.close();
   Gofr.close();
   Gave.close();
 
}




void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}




double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}




void Measure(){ //Properties measurement
  double v, P, t, vij, Pij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp,Pres;

  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;


  Epot.open("Instant/output_epot_"+nomefile+".dat",ios::app);
  Ekin.open("Instant/output_ekin_"+nomefile+".dat",ios::app);
  Temp.open("Instant/output_temp_"+nomefile+".dat",ios::app);
  Etot.open("Instant/output_etot_"+nomefile+".dat",ios::app);
  Pres.open("Instant/output_pres_"+nomefile+".dat",ios::app);
  v = 0.0; //reset observables
  t = 0.0;
  P = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     //update of the histogram of g(r)
     Fill_Histo(dr);


     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       Pij = 48./3./pow(dr,12) - 24.0/3./pow(dr,6);
//Potential energy
       v += vij;
       P += Pij;	
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    walker[iv] = v/(double)npart; //Potential energy
    walker[ik] = t/(double)npart; //Kinetic energy
    walker[it] = (2.0 / 3.0) * t/(double)npart; //Temperature
    walker[ie] = (t+v)/(double)npart; //Total enery
    walker[pr] = rho * walker[it] +  P/vol;


    Epot << walker[iv]  << endl;
    Ekin << walker[ik]  << endl;
    Temp << walker[it] << endl;
    Etot << walker[ie] << endl;
    Pres << walker[pr] << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();

}



	


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf, WriteOld;

  cout << "Print final configuration to file config.final and old.final" << endl << endl;
  WriteConf.open("Conf_final/config_"+nomefile+".final");
  WriteOld.open("Conf_final/old_"+nomefile+".final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteOld << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close(); 
  WriteOld.close();

}





void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
  return;
}




double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}


void Fill_Histo(double dr){
   int bin;
   double dr_box=2.*dr/box*100.;
   bin=floor(dr_box);
   walker[bin+igofr]=walker[bin+igofr]+2.;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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


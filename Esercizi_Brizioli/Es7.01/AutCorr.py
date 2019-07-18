import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os
import pickle
from scipy import optimize

def Correlation (x, delta):
  N=len(x)-delta
  cor=np.zeros(delta)
  for j in range (Delta):
    sum_c=0.
    sum_ave=0.
    sum_1=0.
    sum_2=0.
    for i in range (N):
      sum_c +=x[i]*x[i+j]  
      sum_ave+=x[i+j]
      sum_1 +=x[i]
      sum_2 +=x[i]**2
    sum_1 /=N
    sum_2 /=N
    sum_c /=N
    sum_ave /=N
    cor[j]=(sum_c-(sum_1*sum_ave))/(sum_2-(sum_1)**2)
  return cor


def exponential(x, a, b, c):
  return a * np.exp(-b * x) + c




if not os.path.exists("Data/ACorrel_solid.pkl"):
  Step, Ene_s, Pres_s = np.loadtxt("Data/Inst_solid.0", usecols=(0,1,2), delimiter='\t', unpack='true')
  Delta=1000
  corr_u_s=Correlation(Ene_s,Delta)
  corr_p_s=Correlation(Pres_s,Delta)
  n_step=np.arange(Delta)
  
  output = open('Data/ACorrel_solid.pkl', 'wb')
  pickle.dump((Step,Ene_s,Pres_s,corr_u_s,corr_p_s,n_step), output)
  output.close()

if not os.path.exists("Data/ACorrel_liquid.pkl"):
  Step, Ene_l, Pres_l = np.loadtxt("Data/Inst_liquid.0", usecols=(0,1,2), delimiter='\t', unpack='true')
  Delta=1000
  corr_u_l=Correlation(Ene_l,Delta)
  corr_p_l=Correlation(Pres_l,Delta)
  n_step=np.arange(Delta)
  
  output = open('Data/ACorrel_liquid.pkl', 'wb')
  pickle.dump((Step,Ene_l,Pres_l,corr_u_l,corr_p_l,n_step), output)
  output.close()

if not os.path.exists("Data/ACorrel_gas.pkl"):
  Step, Ene_g, Pres_g = np.loadtxt("Data/Inst_gas.0", usecols=(0,1,2), delimiter='\t', unpack='true')
  Delta=1000
  corr_u_g=Correlation(Ene_g,Delta)
  corr_p_g=Correlation(Pres_g,Delta)
  n_step=np.arange(Delta)
  
  output = open('Data/ACorrel_gas.pkl', 'wb')
  pickle.dump((Step,Ene_g,Pres_g,corr_u_g,corr_p_g,n_step), output)
  output.close()




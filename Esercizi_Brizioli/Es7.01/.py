import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os
import pickle
from scipy import optimize

def error(AV,AV2,n):  # Function for statistical uncertainty estimation
    if n==0:
        return 0
    else:
        return math.sqrt((AV2 - AV**2)/n)


def Media_blk(Step,x)

n_throws=len(Step)


sum_prog=0
su2_prog=0

L=np.arange(0,5100,100,dtype=int)
L[0]=10

n_blk=np.zeros(len(L),dtype=int)
x_blk=np.zeros(len(L))
x_err=np.zeros(len(L))

for t in range(len(L)):
  n_blk[t]=int(n_throws/L[t])
  sum_prog=0
  su2_prog=0
  for i in range(n_blk[t]):
    sum = 0
    for j in range (L[t]):
      k = j+i*L[t]
      sum += x[k]
    sum_prog +=sum/L[t]
    su2_prog +=(sum/L[t])**2
    
  sum_prog/=(n_blk[t])
  su2_prog/=(n_blk[t])

  x_blk[t]=sum_prog
  x_err[t]=error(sum_prog,su2_prog,i)
	
return L, x_blk, x_err


plt.figure()
f, axarr = plt.subplots(1,2)

axarr[0].plot(L,Ene_blk)
axarr[0].set_xlabel('#L')



axarr[1].plot(L,Ene_err)
axarr[1].set_xlabel('#blocks')





plt.figure()
plt.plot(L,Pres_err)


plt.show()



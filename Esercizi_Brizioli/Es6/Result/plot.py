import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os

directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)

points=100
T = np.linspace(0.2,3.0,num=points)
beta = 1/T
J = 1.0
Ns = 50
th = np.tanh(J/T)
thN= th**Ns
ch = 1/th
e = -J*( th + ch*thN )/( 1 + thN )


Temp,u, Err_u = np.loadtxt("Metro_ave.ene.dat", usecols=(0,2,3), delimiter='\t', unpack='true')

Temp1,u1, Err_u1 = np.loadtxt("Gibbs_ave.ene.dat", usecols=(0,2,3), delimiter='\t', unpack='true')


plt.figure()

plt.errorbar(Temp,u,yerr=Err_u,)
plt.errorbar(Temp1,u1,yerr=Err_u1,)
plt.plot(T, e,'r')

plt.title('Ising 1D, internal energy')
plt.xlabel('T')
plt.ylabel('U/N')
plt.show()



import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os

directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)


sigma=0.34*(10**-9)
e_kb=120.
eps=120.*1.380649*(10**-23)
m=39.948*1.66054*(10**-27)
delta=0.0005

Epot = np.loadtxt("Instant/output_epot_solid.dat", usecols=(0), unpack='true')
Ekin = np.loadtxt("Instant/output_ekin_solid.dat", usecols=(0), unpack='true')
Etot = np.loadtxt("Instant/output_etot_solid.dat", usecols=(0), unpack='true')
Temp = np.loadtxt("Instant/output_temp_solid.dat", usecols=(0), unpack='true')
Pres = np.loadtxt("Instant/output_pres_solid.dat", usecols=(0), unpack='true')
x=x= np.arange(len(Epot))


plt.figure()
plt.plot(x,Epot,'r')
plt.plot(x,Ekin,'g')
plt.plot(x,Etot,'b')
plt.plot(x,Temp,'y')

'''
Epot=Epot*eps
Ekin=Ekin*eps
Etot=Etot*eps
Temp=Temp*e_kb
Pres=Pres*eps/(sigma**3)

x=x*delta




Epot_mean = np.loadtxt("ave_epot.out", usecols=(1), unpack='true')
Ekin_mean = np.loadtxt("ave_ekin.out", usecols=(1), unpack='true')
Etot_mean = np.loadtxt("ave_etot.out", usecols=(1), unpack='true')
Temp_mean = np.loadtxt("ave_temp.out", usecols=(1), unpack='true')
Pres_mean = np.loadtxt("ave_pres.out", usecols=(1), unpack='true')
x1= np.arange(len(Epot_mean))*50

plt.plot(x1,Epot_mean,'r')
plt.plot(x1,Ekin_mean,'g')
plt.plot(x1,Etot_mean,'b')
plt.plot(x1,Temp_mean,'y')

Epot_mean=Epot_mean*eps
Ekin_mean=Ekin_mean*eps
Etot_mean=Etot_mean*eps
Temp_mean=Temp_mean*e_kb
Pres_mean=Pres_mean*eps/(sigma**3)
x1= x1*delta

'''

plt.grid(True)


plt.savefig("figure/Lattice_Random_Walk.pdf")


plt.show()

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os

directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)


x, mean, dev_std = np.loadtxt("Lattice_random_walk.txt", usecols=(0,1,2), delimiter='\t', unpack='true')

plt.figure()
#plt.loglog(x,mean)
plt.errorbar(x,mean,yerr=dev_std)
plt.xlabel('Step')
plt.ylabel(r'$\sqrt{\langle r_{N}^{2} \rangle}$')
plt.title('3D Lattice')

plt.grid(True)


plt.savefig("figure/Lattice_Random_Walk.pdf")


x, mean, dev_std = np.loadtxt("Uniform_random_walk.txt", usecols=(0,1,2), delimiter='\t', unpack='true')

plt.figure()
#plt.loglog(x,mean)
plt.errorbar(x,mean,yerr=dev_std)
plt.xlabel('Step')
plt.ylabel(r'$\sqrt{\langle r_{N}^{2} \rangle}$')
plt.title('Continuum')
plt.grid(True)

plt.savefig("figure/Uniform_random_Walk.pdf")
plt.show()

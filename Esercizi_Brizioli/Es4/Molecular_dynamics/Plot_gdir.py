import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os

directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)



bin_solid, mean_solid, error_solid = np.loadtxt("Gofr/Gave_solid.out", usecols=(0,1,2), delimiter='\t', unpack='true')

bin_liquid, mean_liquid, error_liquid = np.loadtxt("Gofr/Gave_liquid.out", usecols=(0,1,2), delimiter='\t', unpack='true')

bin_gas, mean_gas, error_gas = np.loadtxt("Gofr/Gave_gas.out", usecols=(0,1,2), delimiter='\t', unpack='true')

plt.figure()

plt.errorbar(bin_solid,mean_solid,yerr=error_solid)

plt.figure()

plt.errorbar(bin_liquid,mean_liquid,yerr=error_liquid)

plt.figure()

plt.errorbar(bin_gas,mean_gas,yerr=error_gas)

plt.show()

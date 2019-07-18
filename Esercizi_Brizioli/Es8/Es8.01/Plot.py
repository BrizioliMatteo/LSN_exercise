import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os
def Vpot(x):
    return (x**2 - 2.5)*x**2


directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)



index, x1, y1 = np.loadtxt("Prob_distr.out", usecols=(0,1,2), delimiter='\t', unpack='true')

plt.figure()
plt.bar(x1, y1, width=0.05)
plt.show()

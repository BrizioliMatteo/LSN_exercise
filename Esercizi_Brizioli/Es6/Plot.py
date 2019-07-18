import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os

directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)



m = np.loadtxt("Instant_Conf/Metro_instant_T50_h0_Equilibration.0", usecols=(1), delimiter='\t', unpack='true')
Step= np.arange(len(m))

plt.figure()

plt.plot(Step,m)


plt.show()

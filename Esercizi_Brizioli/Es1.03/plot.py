import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os

directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)


x, mean, dev_std = np.loadtxt("PI.txt", usecols=(0,1,2), delimiter='\t', unpack='true')

plt.figure()
plt.errorbar(x,mean,yerr=dev_std)
plt.xlabel('#blocks')
plt.ylabel('$\pi$')
plt.title('Average')
plt.grid(True)

plt.plot([1,100],[3.1415,3.1415])

plt.savefig("figure/PI.pdf")
plt.show()

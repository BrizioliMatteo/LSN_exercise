import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os

directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)


x, mean, dev_std = np.loadtxt("integrale.txt", usecols=(0,1,2), delimiter='\t', unpack='true')
x, mean_p, dev_std_p = np.loadtxt("integrale.txt", usecols=(0,3,4), delimiter='\t', unpack='true')

f, axarr = plt.subplots(1,2, sharey=True)

axarr[0].errorbar(x,mean,yerr=dev_std)
axarr[0].set_title(r'$p(x)=1$')
axarr[0].set_xlabel('#blocks')
axarr[0].set_ylabel('Int')
axarr[0].grid(True)
axarr[0].hlines(y=1.,xmin=1, xmax=100, color='r', linestyle='-')

axarr[1].errorbar(x,mean_p,yerr=dev_std_p)
axarr[1].set_title(r'$p(x)=2(1-x^2)$')
axarr[1].set_xlabel('#blocks')
axarr[1].set_ylabel('Int')
axarr[1].grid(True)
axarr[1].hlines(y=1.,xmin=1, xmax=100,color='r', linestyle='-')

f.suptitle('Integral evaluation')

plt.savefig("figure/Integrale.pdf")
plt.show()

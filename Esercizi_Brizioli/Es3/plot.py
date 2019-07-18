import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os

directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)


x, C_0, C_dev_std = np.loadtxt("Direct.txt", usecols=(0,1,2), delimiter='\t', unpack='true')
x, C_0_d, C_dev_std_d = np.loadtxt("Discret.txt", usecols=(0,1,2), delimiter='\t', unpack='true')


f, axarr = plt.subplots(1,2, sharey=True)

axarr[0].errorbar(x,C_0,yerr=C_dev_std)
axarr[0].set_title('Direct')
axarr[0].set_xlabel('#blocks')
axarr[0].set_ylabel('C[S(0),0]')
axarr[0].grid(True)
axarr[0].hlines(y=14.975790778311286,xmin=1, xmax=100, color='r', linestyle='-')

axarr[1].errorbar(x,C_0_d,yerr=C_dev_std_d)
axarr[1].set_title('Discretized')
axarr[1].set_xlabel('#blocks')
axarr[1].set_ylabel('C[S(0),0]')
axarr[1].grid(True)
axarr[1].hlines(y=14.975790778311286,xmin=1, xmax=100,color='r', linestyle='-')

f.suptitle('Call Options')

plt.savefig("figure/Call_options.pdf")





x, P_0, P_dev_std = np.loadtxt("Direct.txt", usecols=(0,3,4), delimiter='\t', unpack='true')
x, P_0_d, P_dev_std_d = np.loadtxt("Discret.txt", usecols=(0,3,4), delimiter='\t', unpack='true')



f1, axarr1 = plt.subplots(1,2, sharey=True)

axarr1[0].errorbar(x,P_0,yerr=P_dev_std)
axarr1[0].set_title(r'Direct')
axarr1[0].set_xlabel('#blocks')
axarr1[0].set_ylabel('P[S(0),0]')
axarr1[0].grid(True)
axarr1[0].hlines(y=5.4595325819072364,xmin=1, xmax=100, color='r', linestyle='-')

axarr1[1].errorbar(x,P_0_d,yerr=P_dev_std_d)
axarr1[1].set_title(r'Discretized')
axarr1[1].set_xlabel('#blocks')
axarr1[1].set_ylabel('P[S(0),0]')
axarr1[1].grid(True)
axarr1[1].hlines(y=5.4595325819072364,xmin=1, xmax=100,color='r', linestyle='-')

f1.suptitle('Put Options')

plt.savefig("figure/Put_options.pdf")

plt.show()

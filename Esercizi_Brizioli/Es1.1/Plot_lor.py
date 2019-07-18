import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os

directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)

x1,y1= np.loadtxt("Histo_1.txt", usecols=(2,3), delimiter='\t', unpack='true')
x2,y2= np.loadtxt("Histo_2.txt", usecols=(2,3), delimiter='\t', unpack='true')
x3,y3= np.loadtxt("Histo_10.txt", usecols=(2,3), delimiter='\t', unpack='true')
x4,y4= np.loadtxt("Histo_100.txt", usecols=(2,3), delimiter='\t', unpack='true')

f, axarr = plt.subplots(2,2)
f.subplots_adjust(hspace=0.4, wspace=0.4)

axarr[0, 0].bar(x1, y1, width=0.2)
axarr[0, 0].set_title(r'$N=1$')
axarr[0, 1].bar(x2, y2, width=0.2)
axarr[0, 1].set_title(r'$N=2$')
axarr[1, 0].bar(x3, y3, width=0.2)
axarr[1, 0].set_title(r'$N=10$')
axarr[1, 1].bar(x4, y4, width=0.2)
axarr[1, 1].set_title(r'$N=100$')

f.suptitle('Lorentzian Dice')

plt.savefig("figure/Lor.pdf")
plt.show()

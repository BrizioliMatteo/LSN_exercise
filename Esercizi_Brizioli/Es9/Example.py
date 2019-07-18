import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math
import os
import random
import matplotlib.lines as lines

directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)

i, x, y= np.loadtxt("Maps", usecols=(0,1,2), delimiter='\t', unpack='true')


plt.figure()
plt.plot(x, y, linestyle='-',color='k', marker='o',markerfacecolor='blue')
plt.plot(x[0],y[0],linestyle='-',color='k', marker='o',markerfacecolor='red')
#plt.axes().set_aspect('equal')
plt.xlabel('x')
plt.ylabel('y')

plt.grid(True)
plt.show()

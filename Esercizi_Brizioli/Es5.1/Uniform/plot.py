import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

import os

directory= 'figure' # detect the current working directory and print it
if not os.path.exists(directory):
    os.mkdir(directory)

'''
X ,Y, Z = np.loadtxt("output.XYZ_1s.dat", usecols=(0,1,2),  delimiter='\t', unpack='true')
  
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(X, Y, Z, c=Z, marker='.')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.view_init(10, 30)
plt.show()
'''


n ,mean, dev_std = np.loadtxt("output.distr_2p.dat", usecols=(0,1,2),  delimiter='\t', unpack='true')
plt.figure()
plt.errorbar(n,mean,yerr=dev_std)
plt.show()

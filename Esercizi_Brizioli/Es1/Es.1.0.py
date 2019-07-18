import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import math



x, mean, dev_std = np.loadtxt("average.out", usecols=(0,1,2), delimiter='\t', unpack='true')
plt.figure()
plt.errorbar(x,mean-1/2,yerr=dev_std)
plt.xlabel('#throws')
plt.ylabel(r'$\bar{r}-\frac{1}{2}$')
plt.title('Average')
plt.grid(True)

plt.figure()
x, sigma, sigma_dev_std = np.loadtxt("sigma_squared.out", usecols=(0,1,2), delimiter='\t', unpack='true')
plt.errorbar(x,sigma,yerr=sigma_dev_std)
plt.xlabel('#throws')
plt.ylabel(r'$\sigma^{2}$')
plt.title(r'$\sigma^{2}$')
plt.grid(True)

plt.figure()
x, chi = np.loadtxt("chi_squared.out", usecols=(0,1), delimiter='\t', unpack='true')
atteso=np.ones(100)*100
plt.plot(x,chi,'o')
plt.plot(x,atteso,'r',linewidth=2.0)
plt.xlabel('Interval')
plt.ylabel(r'$\chi^{2}$')
plt.title(r'$\chi^{2}$')


plt.show()


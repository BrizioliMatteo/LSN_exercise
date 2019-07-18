import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
import subprocess
from shutil import *
from glob import glob


mu=np.arange(0.7,0.9,0.005)
sigma=np.arange(0.5,0.7,0.005)

with open('input.dat', 'r') as file:
  data = file.readlines()

Ene=10;

for i in range(len(mu)):
  data[0]=str(mu[i])+"\n"
  for j in range(len(sigma)):
    data[1]=str(sigma[j])+"\n"
    with open('input.dat', 'w') as file:
      file.writelines( data )
    cmd= "./Monte_Carlo.exe"
    value = subprocess.call(cmd, shell = True)

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
import subprocess
from shutil import *
from glob import glob


Tem=np.arange(0.5,2.1,0.1)

with open('input.dat', 'r') as file:
  data = file.readlines()

for number in Tem:
  print(number)
  data[0]=str(number)+"\n"
  with open('input.dat', 'w') as file:
    file.writelines( data )
  cmd= "./Monte_Carlo_ISING_1D.exe"
  value = subprocess.call(cmd, shell = True)

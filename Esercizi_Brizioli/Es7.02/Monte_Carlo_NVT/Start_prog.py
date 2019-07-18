import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
import subprocess
from shutil import *
from glob import glob


Tem=[0.8,1.1,1.2]
rho=[1.1,0.8,0.05]
cut_off=[2.2,2.5,5.]
delta=[0.12,0.2,5.]

with open('input.dat', 'r') as file:
  data = file.readlines()

for i in range(3):
  data[0]=str(Tem[i])+"\n"
  data[2]=str(rho[i])+"\n"
  data[3]=str(cut_off[i])+"\n"
  data[4]=str(delta[i])+"\n"
  for l in range (2):
    data[5]=str(10*(l+1))+"\n"
    data[7]=str(l)+"\n"
    with open('input.dat', 'w') as file:
      file.writelines( data )
    cmd= "./Monte_Carlo_NVT.exe"
    value = subprocess.call(cmd, shell = True)

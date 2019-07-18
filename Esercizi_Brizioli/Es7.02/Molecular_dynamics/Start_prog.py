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

with open('input.dat', 'r') as file:
  data = file.readlines()

for i in range(3):
  data[0]=str(Tem[i])+"\n"
  data[2]=str(rho[i])+"\n"
  data[3]=str(cut_off[i])+"\n"
  data[8]=str(0)+"\n"
  with open('input.dat', 'w') as file:
    file.writelines( data )
  cmd= "./MolDyn_NVE.exe"
  value = subprocess.call(cmd, shell = True)

  data[8]=str(1)+"\n"
  for k in range(4):
    with open('input.dat', 'w') as file:
      file.writelines( data )
    cmd= "./MolDyn_NVE.exe"
    value = subprocess.call(cmd, shell = True)

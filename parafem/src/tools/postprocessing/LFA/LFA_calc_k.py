#!/usr/bin/python

# Add path below to .bash_profile (or similar)
#PYTHONPATH=~/PATH_TO_PARAFEM/parafem/bin

import numpy as np
import math
import sys
import os

# Stop .pyc file being created
sys.dont_write_bytecode=True

# Read in file with LFA sample parameters
# Expects fname.py in cwd
cwd = os.getcwd()
fname=sys.argv[1]
sys.path.append(cwd)
module = __import__(fname)

#Define sample variables
l_mm=module.l_mm
r_mm=module.r_mm
m=module.m
Cp=module.Cp
Tmax=module.Tmax

fname=sys.argv[1]+'.ttr2'

#Open file and read data into array
ttr2in = np.fromfile(fname,dtype=float,sep=" ")
length_ttr2in = int(len(ttr2in))
ttr2=np.zeros(shape=(length_ttr2in/2,2),dtype=float)
i=0
for line in range(0,length_ttr2in,2):
    ttr2[i,0]=ttr2in[line]
    ttr2[i,1]=ttr2in[line+1]
    i=i+1

#Find T at half of Tequilibrium
length_ttr2 = int(len(ttr2))
if (Tmax==0.0):
    Tmax=ttr2[length_ttr2-1,1]
T_h=ttr2[0,1]+((Tmax)-ttr2[0,1])/2
#Interpolate to get half-rise time
t_h=np.interp(T_h,ttr2[:,1],ttr2[:,0])

'''
length_ttr2 = int(len(ttr2))
#for row in range(0,length_ttr2):
T=0
i=0

while True:
    T=ttr2[i,1]
    print(T)
    print(i)
    i=i+1
    if (T>=T_h):
        break
'''
#(m^2/s)
l=l_mm*10**-3
alpha=(0.1388*(l**2))/t_h
rho=m/(math.pi*l*(r_mm*10**-3)**2)
k=alpha*Cp*rho

print('alpha = ',alpha,'m^2/s')
print('k = ',k,'W/m K')

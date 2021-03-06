from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os, sys
from collections import Counter
import time
sys.path.append("/home/proj/msc/14/ugrajs/PLUTO41_old/Tools/pyPLUTO")
import pyPLUTO as pp

#Compute how long the simulation takes
start_time = time.time()
#n is an array that stores the size of the simulation domain
n=np.array([256,256,256])

num_bins=1000 #No. of bins of temperature

#constants
UNIT_VELOCITY=(1.e8)
CONST_mp=      1.67262171e-24
CONST_kB=      1.3806505e-16
CONST_mu=      0.5

#Declare all parameters and filenames, file location

filedir="/mnt/lustre/ug4/ugrajs/cooling/higher_k/256/k8-10/"

temp=np.zeros((256,256,256))

for filenumber in xrange(40,41):

    #Load data files using pp.pload

    D=pp.pload(filenumber,w_dir=filedir)

    
    print("--- %s seconds ---" % (time.time() - start_time))

    temp=D.prs/D.rho*(UNIT_VELOCITY*UNIT_VELOCITY)*(CONST_mp*CONST_mu/CONST_kB)    
    temp=np.ndarray.flatten(temp)
    temp=np.sort(temp)
    
    #Make num_bin logarithmic bins of temperature
    temp_binned=np.zeros((num_bins,2))
    Tmax=1e9
    Tmin=1e5
    for i in xrange(0,temp.size):
	if (temp[i]>Tmin):
	    break
    temp=temp[i:]
    index=0
    for i in xrange(0,num_bins):
	temp_binned[i][0]=Tmin*(Tmax/Tmin)**(float(i)/float(num_bins)) 
        for j in xrange(index,temp.size):
	    if(temp[j]>Tmax):
		break
	    elif(temp[j]<temp_binned[i][0]):
		temp_binned[i][1]+=1
	    else:
		index=j+1
		break
    print("--- %s seconds ---" % (time.time() - start_time))
   
    fig, ax = plt.subplots(1)
    fig.set_size_inches(6, 5)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.plot(temp_binned[:,0],temp_binned[:,1])
    plt.xlabel('Temperature (K)')
    plt.ylabel('Number density')
    plt.xlim(1e5,1e8)
    plt.title('PDF as a function of T')
    plt.savefig('256_test.png',dpi=250)
    
 


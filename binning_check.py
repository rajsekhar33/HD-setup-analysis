from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os, sys
from collections import Counter
import time
from scipy.optimize import curve_fit
sys.path.append("/home/rajsekhar/PLUTO41_old/Tools/pyPLUTO")
import pyPLUTO as pp

#Compute how long the simulation takes
start_time = time.time()
#bin_size denotes the division of k space into bins of size 2*pi*bin_size
bin_size=0.75

#n is an array that stores the size of the simulation domain
n=np.array([128,128,128])

#Write k_sq as a functional 3d array, with value at each element given by (nx/2-i)^2+(ny/2-j)^2+(k)^2
k_sq=np.fromfunction(lambda i,j,k:(i-n[0]/2)**2+(j-n[1]/2)**2+k**2,n)

E_k=np.ones(n,dtype=int)
K=np.transpose(np.vstack((np.ndarray.flatten(k_sq),np.ndarray.flatten(E_k))))

#K[][0] stores k, K[][1] stores E(k). 

print("--- %s seconds ---" % (time.time() - start_time))

#The following code adds up E(k) values for the same k^2

c=Counter()
for k, v in K:
#c[k] stores all E(k) values in corresponding k^2 index of the array
	c[k] += v

#K_rad takes non-empty elements of c, and stores the corresponding k(taking square root of k_sq) and E(k) values
K_rad=[[2*np.pi*np.sqrt(var),c[var]] for var in c if var]
K_rad=np.array(K_rad)

#Take bins of a certain size, and add up values corresponding to thse bins

index=0
K_avg=np.empty([0,2],dtype=float)
for i in xrange(2,int((1/bin_size)*np.sqrt(np.sum(np.square(n/2))))):
	k_min=2*bin_size*np.pi*i
	k_max=2*bin_size*np.pi*(i+1)
	E=0
	for j in xrange(index,np.size(K_rad[:,0])):
	    if K_rad[:,0][j]<k_max:
		E+=K_rad[:,1][j]
	        
	    else:
		K_avg=np.append(K_avg,[[k_min,E]],axis=0)
		if(E!=0):
		    index=j
		break
np.savetxt('points_per_bin'+str(bin_size*10)+'_'+str(n[0])+'.txt',K_avg)

#Plot the data 

plt.figure()
plt.plot(K_avg[:,0],K_avg[:,1],'o-')
plt.xlim(0,100)
plt.ylim(0,500)
plt.xlabel('k')
plt.ylabel('Counts')

plt.title('Counts vs k bin size = '+str(bin_size)+', size =' +str(n[0])+'*'+str(n[1])+'*'+str(n[2]))
plt.savefig('counts'+str(bin_size*10)+'_'+str(n[0])+'.png')

print("--- %s seconds ---" % (time.time() - start_time))

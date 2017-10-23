from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os, sys
from collections import Counter
import time
from scipy.optimize import curve_fit
sys.path.append("/home/rajsekhar/PLUTO/Tools/pyPLUTO")
import pyPLUTO as pp

#Compute how long the simulation takes
start_time = time.time()

#average is a quantity that denotes the number of k^2 data points we have averaged over to get Vk_1 and k at the desired point

average=20
#n is an array that stores the size of the simulation domain

n=np.array([128,128,128])

#V1=np.fromfunction(lambda i,j,k:(1/20)*np.cos(np.random.rand()+2*np.pi*20*(-1/2+i/n[0]))+(1/30)*np.cos(np.random.rand()+2*np.pi*30*(-1/2+i/n[0]))+(1/60)*np.cos(np.random.rand()+2*np.pi*60*(-1/2+i/n[0]))+(1/40)*np.cos(np.random.rand()+2*np.pi*40*(-1/2+i/n[0]))+(1/50)*np.cos(np.random.rand()+2*np.pi*50*(-1/2+j/n[0]))+(1/10)*np.cos(np.random.rand()+2*np.pi*(10*(-1/2+k/n[2]))),n)

V1=np.fromfunction(lambda i,j,k:100*np.cos(np.random.rand()+2*np.pi*10*(-1/2+i/n[0]))+400*np.cos(np.random.rand()+2*np.pi*20*(-1/2+i/n[0]))+900*np.cos(np.random.rand()+2*np.pi*30*(-1/2+i/n[0]))+1600*np.cos(np.random.rand()+2*np.pi*40*(-1/2+i/n[0]))+2500*np.cos(np.random.rand()+2*np.pi*50*(-1/2+j/n[0]))+3600*np.cos(np.random.rand()+2*np.pi*(60*(-1/2+k/n[2]))),n)

Vk1=np.fft.fftn(V1,s=n)

Vk_1=np.square(np.abs(np.fft.fftshift(Vk1)))

print("--- %s seconds ---" % (time.time() - start_time))

#Write k_sq as a functional 3d array, with value at each element given by (nx/2-i)^2+(ny/2-j)^2+(nz/2-k)^2
k_sq=np.fromfunction(lambda i,j,k:(-n[0]/2+i)**2+(-n[1]/2+j)**2+(-n[2]/2+k)**2,n)
E_k=np.multiply(Vk_1,k_sq)

#Flatten E_k and k now, and store them in a single  2D array

K=np.transpose(np.vstack((np.ndarray.flatten(k_sq),np.ndarray.flatten(E_k))))

#K[][0] stores k, K[][1] stores Vk_1. 

print("--- %s seconds ---" % (time.time() - start_time))

#The following code adds up Vk_1 values for the same k^2

c=Counter()
for k, v in K:
#c[k] stores all Vk_1 values in corresponding k^2 index of the array
	c[k] += v

#K_rad takes non-empty elements of c, and stores the corresponding k(taking square root of k_sq) and Vk_1 values
K_rad=[[2*np.pi*np.sqrt(var),c[var]] for var in c if var]
K_rad=np.array(K_rad)
#the previous procedure leaves out the value corresponding to k=0, so we insert it
K_rad=np.insert(K_rad,0,c[0],axis=0)

#a denotes the size of the truncated array so that k and Vk_1  can be converted to 2-D arrays for averaging
a=np.shape(K_rad)[0]-np.shape(K_rad)[0]%average

#Truncate the array so that it can be reshaped

K_rad=K_rad[0:a]
K_rad_k=K_rad[:,0].reshape((-1,average))
K_rad_E=K_rad[:,1].reshape((-1,average))



#divide into bins of size average, sum up values inside each bin

K_rad_k=np.mean(K_rad_k,1)
K_rad_E=np.sum(K_rad_E,1)

print("--- %s seconds ---" % (time.time() - start_time))

#Storing the data in a text file
K=np.transpose(np.vstack((K_rad_k,K_rad_E)))
np.savetxt('test_power_spectrum_'+str(n[0])+'.txt',K)

#Plot the data 

fig, ax=plt.subplots()
ax.plot(K_rad_k,K_rad_E)
#ax.plot(K_rad_k,(10**11)*np.power(K_rad_k,4))
plt.xlabel('k')
plt.ylabel('Vk_1')
plt.xlim(0,700)
plt.title('Vk_1 vs k')
plt.savefig('test_spectrum.png')

print("--- %s seconds ---" % (time.time() - start_time))
